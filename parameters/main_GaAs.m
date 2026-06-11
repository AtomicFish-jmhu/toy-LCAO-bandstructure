%% 极简LCAO能带计算模型（本科毕业设计）

% 作者：骥
% 日期：2026.6.9

clear all;
close all;

%% 结构与原子参数输入

primeName = 'GaAs'; 

R_List = [
    0 0 0; 
    0.25 0.25 0.25; 
    ];  

a = 1;  

T = 1 / 0.529 * a * [
    4.62349904e-16 2.87509105e+00 2.87509105e+00;
    2.87509105e+00 0.00000000e+00 2.87509105e+00;
    2.87509105e+00 2.87509105e+00 3.52097106e-16;
];

baseLength = sqrt(T(1, 1)^2 + T(1, 2)^2 + T(1, 3)^2);

% R_L_List = [
%     0, 0, 0;
%     ];

RMax = 2; 
R_L_List = inRange(RMax, T); 

R_List = R_List * T;
R_L_List = R_L_List * T; 

syms x y z real;
r = sqrt(x^2 + y^2 + z^2);
alpha = 0.9;
beta = 1.0;
potential = {
    1, -3 * exp(-alpha * r) / r;
    2, -5 * exp(-beta * r) / r;
};  

fprintf('Calculating stationary atomic orbitals...\n');
phi = {
    1, psai(4, 0, 0, 1);
    1, psai(4, 1, -1, 1);
    1, psai(4, 1, 0, 1);
    1, psai(4, 1, 1, 1);
    2, psai(4, 0, 0, 2);
    2, psai(4, 1, -1, 2);
    2, psai(4, 1, 0, 2);
    2, psai(4, 1, 1, 2);
    }; 

if isempty(gcp('nocreate')) 
    parpool;
end  

%% k路径绘图参数

kpoints = [
    0 0 0;
    0.5, 0. , 0.5;
    0.5 , 0.25, 0.75;
    0.375, 0.375, 0.75 ;
    0 0 0;
    0.5, 0.5, 0.5;
    0.625, 0.25 , 0.625;
    0.5 , 0.25, 0.75;
    0.5, 0.5, 0.5;
    0.375, 0.375, 0.75 ;
];  

k_labels = {'\Gamma', 'X', 'W', 'K', '\Gamma', 'L', 'U', 'W', 'L', 'K'}; 

A = T;  
B = 2*pi*transpose(inv(A)); 
kpoints = kpoints * B; 

div = 1000; 

%% k平面绘图参数

div2 = 300; 

kz = 0; 

%% 生成哈密顿矩阵H与重叠矩阵S的系数张量

tic
[H, S] = HSMatrix(R_List, R_L_List, potential, phi);
toc

name = sprintf('%s_HSMatrix_%s.mat', primeName, datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS'));
save(name, 'H', 'S', 'primeName');

%% 广义特征值求解与能带绘制（一维k路径上）

[X1, Y1, k_indices] = band_structure(H, S, kpoints, div, R_List, R_L_List, phi);

name = sprintf('%s_bands_%s.mat', primeName, datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS'));
save(name, 'X1', 'Y1', 'k_indices', 'k_labels');

fprintf('Plotting diagram...\n');
figure(1)
plot(X1, Y1', '*', 'MarkerSize', 1, 'color', 'b')
set(gca, 'XTick', k_indices);
set(gca, 'XTickLabel', k_labels);

%% 广义特征值求解与能带绘制（二维k平面上）

[X2, Y2, Z2] = band_structure3(H, S, baseLength, div2, kz, R_List, R_L_List, phi);

name = sprintf('%s_bands3_%s.mat', primeName, datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS'));
save(name, 'X2', 'Y2', 'Z2');

fprintf('Plotting diagram...\n');
figure(2);
for k = 1:length(Z2(1, 1, :))
    plot3(X2, Y2, Z2(:, :, k), '*', 'MarkerSize', 0.55)
    hold on;
end
hold off;
axis off;

%% 原子轨道定义

function psi = psai(n, l, m, flag)  
    psi = R(n, l, flag) * Y(l, m);
end

function radial = R(n, l, flag)  
    syms x y z r real;

    if isequal([n, l, flag], [4, 0, 1]) 
        zeta = 1.767;
    elseif isequal([n, l, flag], [4, 1, 1]) 
        zeta = 1.555;
    elseif isequal([n, l, flag], [4, 0, 2])  
        zeta = 2.236;
    elseif isequal([n, l, flag], [4, 1, 2]) 
        zeta = 1.862;
    else
        error('Not included!');
    end
    body = r^(n-1) * exp(-zeta*r);

    body_func = matlabFunction(body, 'Vars', r);
    N_nl = 1 / sqrt(integral(@(r) body_func(r).^2 .* r.^2, 0, Inf));  

    body = subs(body, r, sqrt(x^2 + y^2 + z^2));
    radial = N_nl * body;
end

function angular = Y(l, m) 
    syms x y z real;
    r = sqrt(x^2 + y^2 + z^2);

    index = [l, m];
    switch l
        case 0
            N_lm = 1 / (2*sqrt(pi));
        case 1
            N_lm = (1/2) * sqrt(3/pi);
        case 2
            N_lm = (1/2) * sqrt(15/pi);
    end
    switch true
        case isequal(index, [0, 0])
            body = 1;
        case isequal(index, [1, -1])
            body = y / r;
        case isequal(index, [1, 0])
            body = z / r;
        case isequal(index, [1, 1])
            body = x / r;
        case isequal(index, [2, -2])
            body = x*y / r^2;
        case isequal(index, [2, -1])
            body = y*z / r^2;
        case isequal(index, [2, 0])
            body = (2*z^2 - x^2 - y^2) / (2 * sqrt(3) * r^2);
        case isequal(index, [2, 1])
            body = z*x / r^2;
        case isequal(index, [2, 2])
            body = (x^2 - y^2) / (2 * r^2);  
    end
    angular = N_lm * body;
end