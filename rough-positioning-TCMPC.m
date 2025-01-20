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

qi = AcuRobot.ikine(obj_rough, 'q0', qk', 'tol', 1e-3);

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
itermax = 500;
dUpre = zeros(8 * (NUM_CTRL + 1));
err = 0;
flag = 0;

orien_a = ct.a;
trajectory_3_obj4 = [];
error_3_obj4 = [];
error_3_orientation_obj4 = [];

j1 = 0;
j2 = 0;
% 粗定位
for i = 1:itermax

    if norm(ct.t - O) < 1e-3
        break
    end
    
    P = ct.t;
    error_3_orientation_obj4 = [error_3_orientation_obj4, rad2deg(acos(dot(ct.a , orien_a)/(norm(ct.a)*norm(orien_a))))];
    trajectory_3_obj4 = [trajectory_3_obj4, P];

    norm(TO-P);
    if(norm(TO-P) < 0.5 || flag == 1)
        w = P - TO;
        P_proj = TO + dot(w, Tv) * Tv;
        % PM = P_proj + (O - P_proj) / (10 *NUM_CTRL) * min(j2, 10 * NUM_CTRL);
        % PM = P_proj + (O - P_proj) / (5 *NUM_CTRL);
        if norm(O - P_proj) < 1
            PM = P_proj + (O - P_proj);
        else
            PM = P_proj + (O - P_proj) / (4 *NUM_CTRL);
        end

        Sd = PM - P; % 理想前进方向与距离
        Vv = Pv;
        flag = 1;
        error_3_obj4 = [error_3_obj4, norm(P_proj - P)];
        j2 = j2 + 1;
    elseif(flag == 0)
        w = P - P0;
        P_proj = P0 + dot(w, Pv) * Pv;
        % PM = P_proj + (TO - P_proj) / (10 *NUM_CTRL) * min(j1, 10 * NUM_CTRL);
        % PM = P_proj + (TO - P_proj) / (5 *NUM_CTRL) ;
        if norm(TO - P_proj) < 5
            PM = P_proj + (TO - P_proj);
        else
            PM = P_proj + (TO - P_proj) / (4 *NUM_CTRL);
        end
        Sd = PM - P; % 理想前进方向与距离
        Vv = Tv;
        error_3_obj4 = [error_3_obj4, norm(P_proj - P)];
        j1 = j1 + 1;
    end

    dotvalue = dot(Vv, Sd);
    % fprintf("dotvalue = %f\n", dotvalue);

    [H, g, G, lbA, ubA, psi, omega, theta, Theta, Epsilon, CE, CQ] = mpc2qp_rough_positioning(qk, qi', [], qk_1, uk_1, dUpre, Sd, Vv, dotvalue);
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

    % norm(dt * AcuRobot.jacob0(qk_1', 'trans') * U(1:8, 1) - Sd)
    
    % fprintf("rt - r = %d\n", norm(ct.t - c_T.t) - r);

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
    
    if norm(qk - qi') < 1e-1
        break
    end
    
    % pause(0.1);
end

figure;
x_1 = trajectory_3_obj4(1, :);
y_1 = trajectory_3_obj4(2, :);
z_1 = trajectory_3_obj4(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);
% axis equal;

figure;
plot(0:0.1:i/10-0.1, error_3_obj4, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

figure;
plot(0:0.1:i/10-0.1, error_3_orientation_obj4, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

save('data/rough_expm2.mat', 'trajectory_3_obj4', 'error_3_obj4', 'error_3_orientation_obj4');