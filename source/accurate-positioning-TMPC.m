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

global AcuRobot;
AcuRobot = SerialLink(L,'name','AcuRobot');

qk = [1.63482, 0.0199317, -1.00999, 0.0587652, -2.15094, 0.152961, 0, 0]';

global dt
dt = 0.1;
T = 1;
global NUM_CTRL;
NUM_CTRL = round(T / dt); % number of controls.

curpos = AcuRobot.fkine(qk');

n_seg = NUM_CTRL;
v_segment = divide_line_3d([0, 0, 1], [-0.3, 1.5, 2], n_seg);
v_segment = v_segment(:, 2:end);
% qi_segment = zeros(n_seg, 8);

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
    Q = UnitQuaternion(p);

    trans = Q.SE3;
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

qk_1 = qk;
uk_1 = zeros(8, 1);

Gq_dot = zeros(8 * NUM_CTRL, 8 * (NUM_CTRL + 1));
for i = 1:NUM_CTRL
    for j = 1:i
        Gq_dot(8*(i-1)+1:8*i, 8*(j-1)+1:8*j) = eye(8);
    end
end

ct = AcuRobot.fkine(qk');
dUpre = zeros(8 * (NUM_CTRL + 1));

itermax = 100;
err = 0;


trajectory_2_H_20 = [];
error_2_H_20 = [];
error_orientation_2_H_20 = [];

k = 1;
k_track = 0;
qi = qk';
n_max_track = 30;

% 精定位
for i = 1:itermax
    ctd = ct.double();
    trajectory_2_H_20 = [trajectory_2_H_20, ctd(1:3, 4)];

    error_orientation_2_H_20 = [error_orientation_2_H_20, rad2deg(acos(dot(ct.a, orien_a)))];

    OP = ct.t-c_T.t;
    OPr = r * normalize(OP, 'norm');
    rotate_axis = cross(OPr, dO2);
    anglesum = acos(dot(OPr, dO2)/(r * r));
    h = makehgtform('axisrotate', normalize(rotate_axis, 'norm'), (anglesum/(1 * NUM_CTRL)) * min(i, NUM_CTRL));
    ta = anglesum/(1 * NUM_CTRL) * min(i, NUM_CTRL);
    p = [cos(ta/2), sin(ta/2)*normalize(rotate_axis, 'norm')'];
    Q = UnitQuaternion(p);
    trans2 = Q.SE3;
    OM = h(1:3, 1:3)*OPr;
    Sd = OM - OP; % 理想前进方向与距离

    Vv = OP + 0.5 * Sd;
    
    dotvalue = dot(Vv, Sd);
    % fprintf("dotvalue = %f\n", dotvalue);

    k_track = k_track + 1;
    if k_track >= n_max_track || norm(qk - qi') <= 1e-2
        if k > n_seg
            break;
        end
        qi_1 = qi;
        obj_htm = obj_seg(k);
        try
            qi = AcuRobot.ikine(obj_htm, 'q0', qk', 'tol', 1e-3);
        catch
            qi = AcuRobot.ikine(obj_htm, 'q0', qk', 'tol', 1e-1);
        end
        k = k + 1;
        k_track = 0;
    end

    [H, g, G, lbA, ubA, psi, omega, theta, Theta, Epsilon, CE, CQ] = mpc2qp_accurate_comparative_experiment(qk, qi', d12, qk_1, uk_1, dUpre, Sd, Vv, dotvalue);
    [x,fval,exitflag,iter,lambda,auxOutput] = qpOASES(H, g, G, [], [], lbA, ubA);
    qc = psi * qk + omega * uk_1 + theta * x;
    qcr = reshape(qc, [8, NUM_CTRL]);    

    qcuk = CE * Theta * x - CE * Epsilon + CQ * x;
    qcukr = reshape(qcuk, [16, NUM_CTRL]);
    
    qk_1 = qk;
    qk = qcr(:, 1);
    % qk = qcukr(1:8, 1);
    
    ct = AcuRobot.fkine(qk);
    err = err + abs(norm(ct.t - c_T.t) - r);
    
    U = repmat(uk_1, [NUM_CTRL, 1]) + Gq_dot * x;
    Ur = reshape(U, [8, NUM_CTRL]);

    % norm(dt * AcuRobot.jacob0(qk_1', 'trans') * U(1:8, 1) - Sd)
    
    % fprintf("rt - r = %d\n", norm(ct.t - c_T.t) - r);

    error_2_H_20 = [error_2_H_20, norm(ct.t - c_T.t) - r];

    vttt = dt * AcuRobot.jacob0(qk_1', 'trans') * U(1:8, 1);
    % fprintf("SdTv1 =  %f\n", dot(vttt, Sd));
    % norm(ct.t + dt * AcuRobot.jacob0(qk_1', 'trans') * U(1:8, 1) - c_T.t) - r
    Sd;
    % err = err + abs(norm(ct.t + dt * AcuRobot.jacob0(qk_1', 'trans') * U(1:8, 1) - c_T.t) - r);
    % ct = AcuRobot.fkine(qk);

    epl = Theta * x - Epsilon;
    eplr = reshape(epl, [16, NUM_CTRL]);

    
    uk_1 = U(1:8, 1);
    % uk_1 = Ur(1:8, 1);
    dUpre = x;

    if exitflag ~= 0
        disp("wrong!");
        break;
    end
    
    
    % pause(0.1);
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

x_1 = trajectory_2_H_20(1, :);
y_1 = trajectory_2_H_20(2, :);
z_1 = trajectory_2_H_20(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

figure;
plot(0:0.1:length(error_2_H_20)/10-0.1, error_2_H_20, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

figure;
plot(0:0.1:length(error_orientation_2_H_20)/10-0.1, error_orientation_2_H_20, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

save('data/accurate_expm1.mat', 'trajectory_2_H_20', 'error_2_H_20', 'error_orientation_2_H_20', 'c_T', 'r');


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