clc;
clear;
close all;

% 针灸机器人
L(1) = Link('alpha',pi/2,    'a',0, 'offset',0,    'd',192.5,'standard', 'qlim',[-2*pi,2*pi]);
L(2) = Link('alpha',pi, 'a',266, 'offset',pi/2, 'd',0 ,'standard', 'qlim',[-135, 135] * pi / 180);
L(3) = Link('alpha',pi/2 ,   'a',0,'offset',pi/2,    'd',0, 'standard', 'qlim',[-150, 150] * pi / 180);
L(4) = Link('alpha',-pi/2 ,'a',0, 'offset',0,    'd',324, 'standard', 'qlim',[-2*pi,2*pi]);
L(5) = Link('alpha',pi/2,   'a',0, 'offset',0,    'd',0,'standard', 'qlim',[-147, 147] * pi / 180);
L(6) = Link('alpha',0 ,'a',0, 'offset',0,    'd',155, 'standard', 'qlim',[-2*pi,2*pi]);

L(7) = Link('prismatic', 'a',0, 'offset',45,'alpha',0, 'theta',0, 'qlim',[0,50], 'standard');
L(8) = Link('revolute', 'offset', 0 , 'a', 36.12, 'alpha', 0, 'd', 250, 'qlim',[-2*pi,2*pi], 'standard');

AcuRobot = SerialLink(L,'name','AcuRobot');


% Problem parameters
Np = 10;             % Prediction horizon
n = 6;               % State dimension (joint angles)
m = 8;               % Control input dimension (joint velocities)
dt = 0.1;            % Time step
Q = eye(n+m);          % State cost matrix
Q = blkdiag(eye(3), 1*eye(3), zeros(m));
R = 0.1 * eye(m);    % Input cost matrix
total_time_steps = 100; % Total simulation steps

% Initial conditions
q_current = [1.63482, 0.0199317, -1.00999, 0.0587652, -2.15094, 0.152961, 0, 0]';    % Initial joint angles for Jacobian computation
x_current = AcuRobot.fkine(q_current');  % Initial joint angles

% x_ref_full = repmat(x_current.t - 2, 1, total_time_steps); % Full reference trajectory

n_seg = Np;
v_segment = divide_line_3d([0, 0, 1], [-0.3, 1.5, 2], n_seg);
v_segment = v_segment(:, 2:end);
% qi_segment = zeros(n_seg, 8);

curpos = x_current;
orien_a = zeros(3, 1);
r = 30;  % 针灸长度
t2 = SE3(zeros(3), [0 0 r]);
c_T = curpos * t2;

obj_seg = [];
for i = 1:n_seg
    % 目标位姿计算
    v = normalize(v_segment(:, i)', 'norm');
    z = [0, 0, 1];
    vn = cross(z, v);
    nvn = normalize(vn, 'norm');
    theta = acos(dot(v, z)/(norm(v)*norm(z)));
    p = [cos(theta/2), sin(theta/2)*nvn];
    Q1 = UnitQuaternion(p);

    trans = Q1.SE3;
    trans.t = [0; 0; 0];
    obj = curpos * trans;

    % 构建目标向量


    tl = SE3(zeros(3), -v);
    f_T = curpos * tl;
    dO2 = r * normalize(f_T.t - curpos.t, 'norm');

    dO1 = curpos.t - c_T.t;

    d12 = dO2 - dO1;

    FP = c_T.t + dO2;
    obj.t = FP;

    obj_seg = [obj_seg, obj];

    % qi = AcuRobot.ikine(obj, 'q0', qk', 'tol', 0.081);
    % 
    % qi_segment(i, :) = qi;

    orien_a = obj.a;
end


ki = 1;
x_ref = repmat([obj_seg(n_seg).t; obj_seg(n_seg).torpy()'], Np, 1);

x_ref_init = x_ref;

n_max_track = 30;
k_track = 0;

% x_current = x_current.t;
x_current = [x_current.t; x_current.torpy()'];

u_prev = zeros(m, 1);     % Initial joint velocity


delta_u_opt = zeros(m * Np, 1);

% qpOASES setup
import qpOASES.*;

% Create figure for visualization
figure;
hold on;
% plot3(x_ref_full(1, 1:total_time_steps), x_ref_full(2, 1:total_time_steps), x_ref_full(3, 1:total_time_steps), 'r--', 'LineWidth', 1.5); % Reference trajectory
title('Cartesian Space Trajectory Tracking');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

% Initialize storage for results
x_actual = zeros(n, total_time_steps); % To store actual trajectory


trajectory_comp = [];
error_3_comp = [];
error_3_orientation_comp = [];

x_current1 = x_current;
% Main control loop
for t = 1:total_time_steps  % Number of control iterations
    trajectory_comp = [trajectory_comp, x_current(1:3)];
    ct = AcuRobot.fkine(q_current);
    error_3_orientation_comp = [error_3_orientation_comp, rad2deg(acos(dot(ct.a , orien_a) / (norm(ct.a) * norm(orien_a))))];

    OP = ct.t-c_T.t;
    error_3_comp = [error_3_comp, norm(OP) - r];

    OPr = r * normalize(OP, 'norm');
    rotate_axis = cross(OPr, dO2);
    anglesum = acos(dot(OPr, dO2)/(r * r));
    anglesum = min(1, anglesum);
    h = makehgtform('axisrotate', normalize(rotate_axis, 'norm'), (anglesum/(1 * n_seg) ));
    ta = anglesum/(1 * n_seg) * min(t, n_seg);
    % p = [cos(ta/2), sin(ta/2)*normalize(rotate_axis, 'norm')'];
    % Q2 = UnitQuaternion(p);
    % trans2 = Q2.SE3;
    OM = h(1:3, 1:3)*OPr;
    Sd = OM - OP; % 理想前进方向与距离

    % Store actual trajectory
    x_actual(:, t) = x_current; 
    % Initialize augmented matrices
    n_aug = n + m;                % Augmented state size (state + input)
    A_aug = zeros(n_aug * Np, n_aug);   % Augmented state transition matrix
    B_aug = zeros(n_aug * Np, m * Np);  % Augmented input matrix

    % Temporary variables
    x_temp = x_current;           % Initial state
    q_temp = q_current;           % Initial joint angles
    u_temp = u_prev;              % Initial control input (used for prediction)

    for k = 1:Np
        % Update Jacobian and system matrices at time step k
        J_temp = AcuRobot.jacob0(q_temp', 'rpy');
        % J_temp = J_temp(1:3, :);
        
        % J_temp = J(q_temp);        % Jacobian at current joint angles
        A_k = eye(n);              % Discrete state transition matrix
        B_k = J_temp * dt;         % Input matrix

        % Define the augmented matrices for this step
        A_tilde_k = [A_k, B_k; zeros(m, n), eye(m)]; % Augmented A_k
        B_tilde_k = [B_k; eye(m)];                   % Augmented B_k

        % Update augmented A matrix
        if k == 1
            A_aug((k-1)*n_aug+1:k*n_aug, :) = A_tilde_k;
        else
            A_aug((k-1)*n_aug+1:k*n_aug, :) = A_tilde_k * A_aug((k-2)*n_aug+1:(k-1)*n_aug, :);
        end

        % Update augmented B matrix
        for j = 1:k
            if j == k
                B_aug((k-1)*n_aug+1:k*n_aug, (j-1)*m+1:j*m) = B_tilde_k;
            else
                B_aug((k-1)*n_aug+1:k*n_aug, (j-1)*m+1:j*m) = ...
                    A_tilde_k * B_aug((k-2)*n_aug+1:(k-1)*n_aug, (j-1)*m+1:j*m);
            end
        end

        % Use previous MPC solution to predict joint angles
        if k ~= Np
            u_temp = u_prev + delta_u_opt(k*m+1:(k+1)*m);  % First input update
        else
            % u_temp = u_temp + delta_u_opt((k-1)*m+1:k*m);  % Incrementally update
        end

        % Predict next state
        q_temp = q_temp + u_temp * dt;  % Predict next joint angles

    end

    % Cost function matrices
    H = 2 * (B_aug' * kron(eye(Np), Q) * B_aug + kron(eye(Np), R));  % Quadratic term

    f = 2 * B_aug' * kron(eye(Np), Q) * (A_aug * [x_current; u_prev] - repmat([x_ref(1:6); ones(m, 1)], Np, 1)); % Linear term

    % f = -2 * B_aug' * kron(eye(Np), Q) * ([x_ref; zeros(80, 1)]); % Linear term

    % Matrix construction
    I_m = eye(m);
    I_Np = eye(Np);
    one_Np = ones(Np, 1);

    % Matrix for velocity prediction (Lower triangular matrix)
    T_u = kron(tril(ones(Np)), I_m); % Lower triangular matrix for velocity

    % Matrix for position prediction
    K_q = kron((1:Np)', I_m) * dt; % Scaling matrix for initial velocity

    TTq = zeros(Np, Np);
    for ii = 1:Np
        TTq(ii:Np, ii) = 1:Np-ii+1;
    end
    T_q = kron(TTq, I_m) * dt;

    % T_q = kron(tril(ones(Np)), I_m); % Lower triangular matrix for cumulative position
    % T_q = T_q .* kron((1:Np)', ones(1, Np)) * dt; % Incorporate time step

    u_min = -0.25*[2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 0, 0]';
    u_max = 0.25*[2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 0, 0]';

    q_min = -[2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 0, 0]';
    q_max = [2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 2*pi, 0, 0]';

    lb_u = repmat(u_min, Np, 1) - kron(ones(Np, 1), u_prev);
    ub_u = repmat(u_max, Np, 1) - kron(ones(Np, 1), u_prev);

    lb_q = repmat(q_min, Np, 1) - kron(ones(Np, 1), q_current) - K_q * u_prev;
    ub_q = repmat(q_max, Np, 1) - kron(ones(Np, 1), q_current) - K_q * u_prev;

    % Combine constraints
    G = [T_u; T_q];
    lbG = [lb_u; lb_q];
    ubG = [ub_u; ub_q];

    T_deltau = eye(m * Np, m * Np);
    lb_deltau = repmat(-0.1 * ones(m, 1), Np, 1);
    ub_deltau = repmat(0.1 * ones(m, 1), Np, 1);

    G = [G; T_deltau];
    lbG = [lbG; lb_deltau];
    ubG = [ubG; ub_deltau];

    % % 末端姿态约束
    % JT = AcuRobot.jacob0(q_current');
    % A_eq = JT(4:5, :) * [eye(m), zeros(m, (Np-1)*m)];
    % b_eq = [0, 0]' - JT(4:5, :) * u_prev ;
    % 
    % G = [G; A_eq];
    % lbG = [lbG; b_eq];
    % ubG = [ubG; b_eq];

    % 末端速度约束1*m*Np + 3
    JT = AcuRobot.jacob0(q_current', 'rpy');
    V_eq = JT(1:3, :) * [eye(m), zeros(m, (Np-1)*m)] * dt;
    vb_eq = Sd - JT(1:3, :) * u_prev * dt;

    G = [G; V_eq];
    lbG = [lbG; vb_eq];
    ubG = [ubG; vb_eq];

    % introduce slak factor
    % multi slack
    G_s = [G, eye(3*m*Np + 3); -G, -eye(3*m*Np + 3); zeros(3*m*Np + 3, 80), eye(3*m*Np + 3)];
    lbG_s = [lbG; -ubG; zeros(3*m*Np + 3, 1)];
    ubG_s = [ubG; -lbG; inf(3*m*Np + 3, 1)];

    rou = 10;
    H_s = blkdiag(H, 2*rou*eye(3*m*Np + 3));
    f_s = [f; zeros(3*m*Np + 3, 1)];

    % G_s = [T_u, zeros(1*m*Np, 2*m*Np + 5); 
    %     T_q, zeros(1*m*Np, 2*m*Np + 5); 
    %     [T_deltau; A_eq; V_eq], eye(2*m*Np + 5);
    %     -[T_deltau; A_eq; V_eq], -eye(2*m*Np + 5);
    %     zeros(2*m*Np + 5, 80), eye(2*m*Np + 5)];
    % 
    % lbG_s = [lb_u; lb_q; [lb_deltau; b_eq; vb_eq]; -[ub_deltau; b_eq; vb_eq]; zeros(2*m*Np + 5, 1)];
    % ubG_s = [ub_u; ub_q; [ub_deltau; b_eq; vb_eq]; -[lb_deltau; b_eq; vb_eq]; inf(2*m*Np + 5, 1)];
    % 
    % rou = 10;
    % H_s = blkdiag(H, 2*rou*eye(2*m*Np + 5));
    % f_s = [f; zeros(2*m*Np + 5, 1)];

    % G_s = [T_q, zeros(1*m*Np, 2*m*Np + 5);
    %     [T_u; T_deltau; A_eq; V_eq], eye(2*m*Np + 5);
    %     -[T_u; T_deltau; A_eq; V_eq], -eye(2*m*Np + 5);
    %     zeros(2*m*Np + 5, 80), eye(2*m*Np + 5)];
    % 
    % lbG_s = [lb_q; [lb_u; lb_deltau; b_eq; vb_eq]; -[ub_u; ub_deltau; b_eq; vb_eq]; zeros(2*m*Np + 5, 1)];
    % ubG_s = [ub_q; [ub_u; ub_deltau; b_eq; vb_eq]; -[lb_u; lb_deltau; b_eq; vb_eq]; inf(2*m*Np + 5, 1)];
    % 
    % rou = 10;
    % H_s = blkdiag(H, 2*rou*eye(2*m*Np + 5));
    % f_s = [f; zeros(2*m*Np + 5, 1)];
    
    % % single slack
    % G_s = [G, ones(3*m*Np + 3, 1); -G, -ones(3*m*Np + 3, 1); zeros(1, 80), 1];
    % lbG_s = [lbG; -ubG; 0];
    % ubG_s = [ubG; -lbG; inf];
    % 
    % rou = 0.5;
    % H_s = blkdiag(H, 2*rou*1);
    % f_s = [f; 0];

    % JT = AcuRobot.jacob0(q_current', 'eul');
    % A_eq = JT(1:3, :) * [eye(m), zeros(m, (Np-1)*m)] * dt;
    % b_eq = [-0.01, 0, 0]' - JT(1:3, :) * u_prev * dt;
    % 
    % G = [G; A_eq];
    % lbG = [lbG; b_eq];
    % ubG = [ubG; b_eq];


    % Solve QP problem using qpOASES
    options = qpOASES_options('default');
    options.printLevel = 1;
    % qp_solver = qpOASES();
    
    % Constraints (no inequality for simplicity, add as needed)
    lb = -0.1 * ones(m * Np, 1);  % Lower bounds for delta u
    ub = 0.1 * ones(m * Np, 1);   % Upper bounds for delta u

    % Solve QP
    % [delta_u_opt, obj_val, exit_flag] = qp_solver(H, f, [], [], [], lb, ub);

    % [delta_u_opt, obj_val, exit_flag] = qpOASES(H, f, G, [], [], lbG, ubG, options);

    [delta_u_opt, obj_val, exit_flag] = qpOASES(H_s, f_s, G_s, [], [], lbG_s, ubG_s, options);

    if exit_flag ~= 0
        disp('QP solver failed.');
        break;
    end

    % Apply first control input
    delta_u_1 = delta_u_opt(1:m);  % First input increment
    u_prev = u_prev + delta_u_1;   % Update velocity control input

    % Update Cartesian state
    JT = AcuRobot.jacob0(q_current', 'rpy');
    x_current1 = x_current1 + JT * u_prev * dt;  % Update Cartesian state

    % Update joint angles (intermediate variable for Jacobian)
    q_current = q_current + u_prev * dt;  % Update joint angles

    
    ti = AcuRobot.fkine(q_current);
    x_current = [ti.t; ti.torpy()'];
    
    % Log or display results
    disp(['Iteration ', num2str(t), ', obj_val: ', num2str(obj_val)]);
    disp(['Current state: ', mat2str(x_current')]);

    k_track = k_track + 1;
    % if k_track >= n_max_track || norm(x_current(1:6) - x_ref(1:6)) < 1e-1
    %     ki = ki + 1;
    %     if ki <= n_seg
    %         x_ref = repmat([obj_seg(ki).t; obj_seg(ki).torpy()'], Np, 1);
    %         k_track = 0;
    %     end
    % 
    % end

    if norm(x_current - [obj_seg(n_seg).t; obj_seg(n_seg).torpy()']) < 1e-1
        break
    end
end


center = c_T.t;
radius = r;

% Create a sphere
[X, Y, Z] = sphere(100); % Increase the resolution by using 50 points

% Scale and shift the sphere
X = radius * X + center(1);
Y = radius * Y + center(2);
Z = radius * Z + center(3);

% Plot the sphere
figure;
h = surf(X, Y, Z);

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Transparent Sphere');
axis equal; % Ensures that the scale of the axes are equal

% Set color and transparency
shading interp; % Smooth shading

set(h, 'FaceColor', [0.94118, 0.9755, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

hold on;

x_1 = trajectory_comp(1, :);
y_1 = trajectory_comp(2, :);
z_1 = trajectory_comp(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

figure;
plot(0:0.1:length(error_3_comp)/10-0.1, error_3_comp, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

figure;
plot(0:0.1:length(error_3_orientation_comp)/10-0.1, error_3_orientation_comp, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

save('data/accurate_expm3.mat', 'trajectory_comp', 'error_3_comp', 'error_3_orientation_comp');


function segments = divide_line_3d(p1, p2, n)
    % p1: 第一个点的坐标，列向量 [x1; y1; z1]
    % p2: 第二个点的坐标，列向量 [x2; y2; z2]
    % n: 需要划分的线段数
    % segments: 返回所有线段端点的坐标，大小为 3 x (n+1) 的矩阵

    % 预分配存储线段端点的矩阵
    segments = zeros(3, n+1);
    
    % 计算每个线段在 x, y 和 z 方向上的增量
    dx = (p2(1) - p1(1)) / n;
    dy = (p2(2) - p1(2)) / n;
    dz = (p2(3) - p1(3)) / n;
    
    % 计算每个端点的坐标
    for i = 0:n
        segments(1, i+1) = p1(1) + i * dx;
        segments(2, i+1) = p1(2) + i * dy;
        segments(3, i+1) = p1(3) + i * dz;
    end
end