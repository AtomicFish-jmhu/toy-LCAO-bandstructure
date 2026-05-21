%% 本脚本接收main.m输出的矩阵H和S，计算并绘制沿指定k空间路径的能带

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
];  % 存放高对称k点的分数坐标，在原胞基矢对应的倒格基矢下表达

k_labels = {'\Gamma', 'X', 'W', 'K', '\Gamma', 'L', 'U', 'W', 'L', 'K'};  % 节点标签

A = T;  
% A = 1 / 0.529 * a * [
% 
% ];  % 原胞基矢的直角坐标，一排一个基矢。以防你的T用的是惯用基，重新输入一遍

B = 2*pi*transpose(inv(A));  % 原胞基矢对应的倒格基矢，直角坐标表示，一排一个

kpoints = kpoints * B;  % 将高对称点分数坐标转成直角坐标

div = 1000;  % 单位k空间长度上的划分数

[k_indices, X, Y] = band(H, S, kpoints, div);

currentTime = datetime('now');
timeStr = datestr(currentTime, 'yyyy-mm-dd_HH-MM-SS');
name = sprintf('bands(k line)_%s_%s.mat', filename, timeStr);
save(name, 'X', 'Y', 'k_indices', 'k_labels'); 

plot(X, Y, '*', 'MarkerSize', 1, 'color', 'b')
set(gca, 'XTick', k_indices);
set(gca, 'XTickLabel', k_labels);

function [k_indices, X, Y] = band(H, S, kpoints, div)  % 生成沿指定k空间高对称路径上的能带

% 先由高对称节点坐标生成路径点列

    fprintf('Generating dots sequence...');
    
    N = length(kpoints(:, 1));  % 拐点数
    
    path_length = zeros(N-1, 1);  % 每段长，存入数组
    for i = 1:(N-1)
        path_length(i) = sqrt((kpoints(i, 1) - kpoints(i+1, 1))^2 + ...
                              (kpoints(i, 2) - kpoints(i+1, 2))^2 + ...
                              (kpoints(i, 3) - kpoints(i+1, 3))^2);
    end
    
    points_num = zeros(N-1, 1);  % 每段上的取样数目，不包括最后一段最后一个点
    total_points_num = 0;  % 总取样数目
    for i = 1:(N-1)
        points_num(i) = ceil(path_length(i) * div);
        total_points_num = total_points_num + points_num(i);
    end
    total_points_num = total_points_num + 1;  % 补上最后一个点
    
    kpath = zeros(total_points_num, 3);  % 记录所有取样点坐标的数组
    trace = 0;  % 已取样的点总数
    for i = 1:(N-1)
        for j = 1:3
            temp = linspace(kpoints(i, j), kpoints(i+1, j), points_num(i)+1).';
            kpath(trace + 1 : (trace + points_num(i)), j) = temp(1 : points_num(i));
        end
        trace = trace + points_num(i);
    end
    kpath(total_points_num, 1:3) = kpoints(N, 1:3);  % 补上最后一个节点
    
    k_indices = zeros(1, N);  % 存放各个拐点的索引值，用于绘图时标注
    k_indices(1) = 1;
    for i = 2:N
        k_indices(i) = k_indices(i-1) + points_num(i-1);
    end 
    
% 计算指定路径kpath上各点的能量
    
    band_num = length(H(:, 1));
    X = linspace(1, total_points_num, total_points_num);
    Y = zeros(total_points_num, band_num);
    
    syms k_x k_y k_z real;
    
    fprintf('\nSerializing matlab functions...');
    H = matlabFunction(H, 'var', [k_x, k_y, k_z]); 
    S = matlabFunction(S, 'var', [k_x, k_y, k_z]);
        
    fprintf('\nCalculating bands structure...');
    for i = 1:total_points_num
        H_val = H(kpath(i, 1), kpath(i, 2), kpath(i, 3));
        S_val = S(kpath(i, 1), kpath(i, 2), kpath(i, 3));
        solution = eig(H_val, S_val);
        Y(i, :) = real(solution); 
        if mod(i, floor(total_points_num * 0.1)) == 0
            fprintf('■');
        end
    end
    fprintf(' Complete!\n');
end