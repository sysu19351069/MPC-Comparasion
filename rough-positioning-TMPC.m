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

% qk = [1.63482, 0.0199317, -1.00999, 0.0587652, -2.15094, 0.152961, 0, 0]';

% 1、起始点为[29.9135, 266.1760, 60.4643]
% qk = [-4.5221, 0.2793, 1.6748, -0.1403, 1.7914, -2.6356, 0, 0]';

% 2、起始点为[212.1069, -383.3628, 101.9846]
% qk = [-1.2637, -1.2590, -0.6789, -2.9400, -2.5430, 4.7357, 0, 0]';

% 3、起始点为[45.1866, 96.6825, 66.6616]
qk = [0.3140   -1.3949   -1.9689    2.9790    2.5205   -1.5637         0         0]';

global dt
dt = 0.1;
T = 1;
global NUM_CTRL;
NUM_CTRL = round(T / dt); % number of controls.

curpos = AcuRobot.fkine(qk');

obj_rough = curpos;
obj_rough.t = [500; 110; -100];

% obj1
% obj_rough.t = [300; 300; -50];

% obj2
% obj_rough.t = [400; -250; -200];

% obj3
% obj_rough.t = [450; 80; -100];

% obj4
% obj_rough.t = [580; 100; -150];

qi = AcuRobot.ikine(obj_rough, 'q0', qk', 'tol', 0.081);

P0 = curpos.t;
O = obj_rough.t;
TO = [O(1:2); P0(3)];
Pv = normalize(TO - P0, 'norm');
Tv = normalize(O - TO, 'norm');

qk_1 = qk;
uk_1 = zeros(8, 1);

Gq_dot = zeros(8 * NUM_CTRL, 8 * (NUM_CTRL + 1));
for i = 1:NUM_CTRL
    for j = 1:i
        Gq_dot(8*(i-1)+1:8*i, 8*(j-1)+1:8*j) = eye(8);
    end
end

ct = AcuRobot.fkine(qk);
itermax = 1500;
dUpre = zeros(8 * (NUM_CTRL + 1));
err = 0;
flag = 0;

orien_a = ct.a;
trajectory_comp = [];
error_3_comp = [];
error_3_orientation_comp = [];

j1 = 0;
j2 = 0;

% dis = 5;
dis = norm(TO - P0) / (4*NUM_CTRL);
n1 = floor(norm(TO - P0)/dis);

dis = norm(O - TO) / (4*NUM_CTRL);
n2 = floor(norm(O - TO)/dis);
segments1 = divide_line_3d(P0, TO, n1);
segments2 = divide_line_3d(TO, O, n2);

segments = [segments1(:, 2:end), segments2(:, 2:end)];
n = n1 + n2;

% segments = O;
% n = 1;

obj_point = P0;
k = 1;

qi = qk';

err_qi = [];

k_track = 0;
n_max_track = 30;

% 粗定位
for i = 1:itermax
    
    P = ct.t;
    error_3_orientation_comp = [error_3_orientation_comp, rad2deg(acos(dot(ct.a , orien_a)/(norm(ct.a)*norm(orien_a))))];
    trajectory_comp = [trajectory_comp, P];

    norm(TO-P);

    w = P - P0;
    P_proj1 = P0 + dot(w, Pv) * Pv;
    
    w = P - TO;
    P_proj2 = TO + dot(w, Tv) * Tv;

    if norm(P_proj1 - P) < norm(P_proj2 - P)
        error_3_comp = [error_3_comp, norm(P_proj1 - P)];
    else
        error_3_comp = [error_3_comp, norm(P_proj2 - P)];
    end

    % if norm(P - obj_point) < 2
    %     if k > n
    %         break;
    %     end
    %     qi_1 = qi;
    %     obj_point = segments(:, k);
    %     obj_htm = ct;
    %     obj_htm.t = obj_point;
    %     qi = AcuRobot.ikine(obj_htm, 'q0', qk', 'tol', 5e-3);
    %     k = k + 1;
    %     err_qi = [err_qi, norm(qi_1 - qi)];
    % end
    
    k_track = k_track + 1;
    if k_track >= n_max_track || norm(qk - qi') <= 1e-3
        if k > n
            break;
        end
        qi_1 = qi;
        obj_point = segments(:, k);
        obj_htm = ct;
        obj_htm.t = obj_point;
        try
            qi = AcuRobot.ikine(obj_htm, 'q0', qk', 'tol', 1e-3);
            err_qi = [err_qi, norm(qi_1 - qi)];
        catch
            qi = AcuRobot.ikine(obj_htm, 'q0', qk', 'tol', 2e-1);
            err_qi = [err_qi, norm(qi_1 - qi)];
        end
        k = k + 1;
        k_track = 0;
    end

    [H, g, G, lbA, ubA, psi, omega, theta, Theta, Epsilon, CE, CQ] = mpc2qp_comparative_experiment(qk, qi', [], qk_1, uk_1, dUpre, [], [], []);
    [x,fval,exitflag,iter,lambda,auxOutput] = qpOASES(H, g, G, [], [], lbA, ubA);
    qc = psi * qk + omega * uk_1 + theta * x;
    qcr = reshape(qc, [8, NUM_CTRL]);    

    qcuk = CE * Theta * x - CE * Epsilon + CQ * x;
    qcukr = reshape(qcuk, [16, NUM_CTRL]);
    
    qk_1 = qk;
    qk = qcr(:, 1);
    % qk = qcukr(1:8, 1);
    
    ct = AcuRobot.fkine(qk);
    
    U = repmat(uk_1, [NUM_CTRL, 1]) + Gq_dot * x;
    Ur = reshape(U, [8, NUM_CTRL]);

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

figure;
x_1 = trajectory_comp(1, :);
y_1 = trajectory_comp(2, :);
z_1 = trajectory_comp(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

% zticks([-100, 0, 100]);

set(gca, 'FontSize', 20,  'FontName', 'Times');

xlabel('x(mm)', 'FontSize', 20, 'FontName', 'Times');
ylabel('y(mm)', 'FontSize', 20, 'FontName', 'Times');
zlabel('z(mm)', 'FontSize', 20, 'FontName', 'Times');

lgd = legend({'start1'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

% axis equal;

figure;
set(gca, 'FontSize', 20,  'FontName', 'Times');
plot(0:0.1:i/10-0.1, error_3_comp, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);
ylabel('Right angle tracking error/mm', 'FontSize', 18, 'FontName', 'Times');

xlabel('t/s', 'FontSize', 20, 'FontName', 'Times');

lgd = legend({'Start1\_tte'}, 'Location', 'best');
lgd.FontSize = 16; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

figure;
plot(0:0.1:i/10-0.1, error_3_orientation_comp, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

save('data/rough_expm1.mat', 'trajectory_comp', 'error_3_comp', 'error_3_orientation_comp');


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