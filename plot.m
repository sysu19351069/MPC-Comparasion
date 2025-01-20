close all;
clear;

%% 直角轨迹跟踪
figure(1);
hold on;

load 'data/rough_expm1.mat';
x_1 = trajectory_comp(1, :);
y_1 = trajectory_comp(2, :);
z_1 = trajectory_comp(3, :);

plot3(x_1, y_1, z_1, '--', 'LineWidth', 1, 'Color', [0, 0, 1]);

% zticks([-100, 0, 100]);

set(gca, 'FontSize', 20,  'FontName', 'Times');

xlabel('x(mm)', 'FontSize', 20, 'FontName', 'Times');
ylabel('y(mm)', 'FontSize', 20, 'FontName', 'Times');
zlabel('z(mm)', 'FontSize', 20, 'FontName', 'Times');

figure(2);
hold on;
set(gca, 'FontSize', 20,  'FontName', 'Times');
plot(0:0.1:length(error_3_comp)/10-0.1, error_3_comp, '-', 'LineWidth', 1, 'Color', [0, 0, 1]);
ylabel('Right angle tracking error/mm', 'FontSize', 18, 'FontName', 'Times');

xlabel('t/s', 'FontSize', 20, 'FontName', 'Times');

figure(3);
hold on;

plot(0:0.1:length(error_3_orientation_comp)/10-0.1, error_3_orientation_comp, '-', 'LineWidth', 1, 'Color', [0, 0, 1]);
set(gca, 'FontSize', 20,  'FontName', 'Times');

xlabel('t/s', 'FontSize', 20, 'FontName', 'Times');
ylabel('Orientation error /^\circ', 'FontSize', 18, 'FontName', 'Times');

%
figure(1);
hold on;

load 'data/rough_expm3.mat';
x_1 = trajectory_comp(1, :);
y_1 = trajectory_comp(2, :);
z_1 = trajectory_comp(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [0, 1, 0]);

figure(2);
hold on;
plot(0:0.1:length(error_3_comp)/10-0.1, error_3_comp, '-', 'LineWidth', 2, 'Color', [0, 1, 0]);

figure(3);
hold on;

plot(0:0.1:length(error_3_orientation_comp)/10-0.1, error_3_orientation_comp, '-', 'LineWidth', 1, 'Color', [0, 1, 0]);

%
figure(1);
hold on;

load 'data/rough_expm2.mat';
x_1 = trajectory_3_obj4(1, :);
y_1 = trajectory_3_obj4(2, :);
z_1 = trajectory_3_obj4(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 2, 'Color', [1, 0, 0]);

lgd = legend({'TMPC', 'LTV-MPC', 'TCMPC'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

figure(2);
hold on;
plot(0:0.1:length(error_3_obj4)/10-0.1, error_3_obj4, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);
lgd = legend({'TMPC', 'LTV-MPC', 'TCMPC'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

figure(3);
hold on;

plot(0:0.1:length(error_3_orientation_obj4)/10-0.1, error_3_orientation_obj4, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);
lgd = legend({'TMPC', 'LTV-MPC', 'TCMPC'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

%% 球面轨迹跟踪
figure(4);
hold on;

load 'data/accurate_expm1.mat';
center = c_T.t;
radius = r;

% Create a sphere
[X, Y, Z] = sphere(100); % Increase the resolution by using 50 points

% Scale and shift the sphere
X = radius * X + center(1);
Y = radius * Y + center(2);
Z = radius * Z + center(3);

% Plot the sphere
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

x_1 = trajectory_2_H_20(1, :);
y_1 = trajectory_2_H_20(2, :);
z_1 = trajectory_2_H_20(3, :);

plot3(x_1, y_1, z_1, '--', 'LineWidth', 1, 'Color', [0, 0, 1]);

% zticks([-100, 0, 100]);

set(gca, 'FontSize', 20,  'FontName', 'Times');

xlabel('x(mm)', 'FontSize', 20, 'FontName', 'Times');
ylabel('y(mm)', 'FontSize', 20, 'FontName', 'Times');
zlabel('z(mm)', 'FontSize', 20, 'FontName', 'Times');

% lgd = legend({'start1'}, 'Location', 'best');
% lgd.FontSize = 18; % 设置字体大小
% lgd.FontName = 'Times'; % 设置字体名称为Times

% axis equal;

figure(5);
hold on;
set(gca, 'FontSize', 20,  'FontName', 'Times');
plot(0:0.1:length(error_2_H_20)/10-0.1, error_2_H_20, '-', 'LineWidth', 1, 'Color', [0, 0, 1]);
ylabel('Sphere tracking error/mm', 'FontSize', 18, 'FontName', 'Times');

xlabel('t/s', 'FontSize', 20, 'FontName', 'Times');

% lgd = legend({'Start1\_tte'}, 'Location', 'best');
% lgd.FontSize = 16; % 设置字体大小
% lgd.FontName = 'Times'; % 设置字体名称为Times

figure(6);
hold on;

plot(0:0.1:length(error_orientation_2_H_20)/10-0.1, error_orientation_2_H_20, '-', 'LineWidth', 1, 'Color', [0, 0, 1]);
set(gca, 'FontSize', 20,  'FontName', 'Times');

xlabel('t/s', 'FontSize', 20, 'FontName', 'Times');
ylabel('Orientation error /^\circ', 'FontSize', 18, 'FontName', 'Times');


%
figure(4);
hold on;

load 'data/accurate_expm3.mat';
x_1 = trajectory_comp(1, :);
y_1 = trajectory_comp(2, :);
z_1 = trajectory_comp(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [0, 1, 0]);

figure(5);
hold on;
plot(0:0.1:length(error_3_comp)/10-0.1, error_3_comp, '-', 'LineWidth', 2, 'Color', [0, 1, 0]);

figure(6);
hold on;

plot(0:0.1:length(error_3_orientation_comp)/10-0.1, error_3_orientation_comp, '-', 'LineWidth', 1, 'Color', [0, 1, 0]);

%
figure(4);
hold on;

load 'data/accurate_expm2.mat';
x_1 = trajectory_2_H_20(1, :);
y_1 = trajectory_2_H_20(2, :);
z_1 = trajectory_2_H_20(3, :);

plot3(x_1, y_1, z_1, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);

lgd = legend({'Sphere', 'TMPC', 'LTV-MPC', 'TCMPC'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

figure(5);
hold on;
plot(0:0.1:length(error_2_H_20)/10-0.1, error_2_H_20, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);
lgd = legend({'TMPC', 'LTV-MPC', 'TCMPC'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times

figure(6);
hold on;

plot(0:0.1:length(error_orientation_2_H_20)/10-0.1, error_orientation_2_H_20, '-', 'LineWidth', 1, 'Color', [1, 0, 0]);
lgd = legend({'TMPC', 'LTV-MPC', 'TCMPC'}, 'Location', 'best');
lgd.FontSize = 18; % 设置字体大小
lgd.FontName = 'Times'; % 设置字体名称为Times
