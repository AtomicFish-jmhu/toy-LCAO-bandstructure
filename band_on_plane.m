%% 本脚本接收main.m输出的矩阵H和S，计算指定k_x k_y平面上的能带，输出定义域X, Y和能量值Z。绘图需与band_on_plane_plot连用

range = 6*pi/baseLength;  % 展示的k空间范围
div = 300;  % k空间取样点数目
kz = 0;  % 默认绘图固定k_z

[X, Y, Z] = band(H, S, filename, range, div, kz); 

function [X, Y, Z] = band(H, S, filename, range, div, kz)  % 生成指定k_x k_y平面上的能带
    syms k_x k_y k_z real;
    
    kx_vals = linspace(-range/2, range/2, div);
    ky_vals = linspace(-range/2, range/2, div);
    [X, Y] = meshgrid(kx_vals, ky_vals);
    
    band_num = length(H(:,1));
    Z = zeros(div, div, band_num);
    
    H = matlabFunction(H, 'var', [k_x, k_y, k_z]);  % 数值函数比符号函数的计算快很多。不要在循环体中用符号函数计算
    S = matlabFunction(S, 'var', [k_x, k_y, k_z]);
    
    fprintf('Calculating bands structure...');
    for i = 1:div
        parfor j = 1:div
            H_val = H(kx_vals(i), ky_vals(j), kz);
            S_val = S(kx_vals(i), ky_vals(j), kz);
            solution = eig(H_val, S_val);  % 求解由方程 Hx = ESx 定义的广义特征值
            Z(i, j, :) = real(solution); 
        end
        if mod(i, div * 0.1) == 0
            fprintf('■');  % 在最初的版本中，函数输出“吱！”报告进度
        end
    end
    fprintf(' Complete!\n');
    
    currentTime = datetime('now');
    timeStr = datestr(currentTime, 'yyyy-mm-dd_HH-MM-SS');
    name = sprintf('bands(k plane)_%s_%s.mat', filename, timeStr);
    save(name, 'X', 'Y', 'Z', 'range', 'div', 'kz');
end