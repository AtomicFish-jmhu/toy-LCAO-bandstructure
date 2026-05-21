%% 基于全积分LCAO方法的能带计算模型（本科毕业设计）
% 作者：骥
% 日期：2026.4
% 提醒：这是 TMatrix 版本！未经充分测试且不再维护。

%% 在此脚本中输入计算参数，包括单原子势场，以及参与计算的原子坐标和原子轨道。运行main.m后输出哈密顿矩阵H和重叠矩阵S。

clear all;
close all;

R_List = [
    0 0 0;  % Ga
    0.25 0.25 0.25;  % As
    ];  % 存放一个原胞中各原子的坐标，用原胞基矢表示

% R_L_List = [
%     0, 0, 0;  % 中心原胞
% ];  % 存放要遍历的格矢坐标，用原胞基矢表示。仅在你想要手动输入时启用

a = 1;  
% 晶格伸缩系数，用于调试（a太小将引发发散困难；a太大将使能级简并，观察不到能带结构）

A = 1 / 0.529 * a * [
    4.62349904e-16 2.87509105e+00 2.87509105e+00;
    2.87509105e+00 0.00000000e+00 2.87509105e+00;
    2.87509105e+00 2.87509105e+00 3.52097106e-16;
];  
% 原胞基矢的直角坐标，横着摆放，一排一个基矢
% 1 / 0.529 为从 Å 转到原子单位制时做的转换（1 a.u. = 1 波尔半径 a_0 ≈ 0.529 Å）

baseLength = sqrt(A(1, 1)^2 + A(1, 2)^2 + A(1, 3)^2);

RMax = 2;  % 截断半径，以第一个原胞基矢长度为单位1
R_L_List = inRange(RMax, A);  % 添加截断半径以内的格矢

R_List = R_List * A;
R_L_List = R_L_List * A;  % 将原子位矢和原胞格矢转为直角坐标

syms x y z real;
r = sqrt(x^2 + y^2 + z^2);
potential = {
    1, -3 / r * exp(-1.0 * r);
    2, -5 / r * exp(-0.5 * r);
    };  
% 单原子势场，可以改成别的球对称势 V(r)
% 前面的数字表示势场属于 R_List 中第几排坐标代表的原子

fprintf('Calculating stationary atomic orbitals...\n');
phi = {
    1, psai(4, 0, 0, 1), -0.424;
    1, psai(4, 1, -1, 1), -0.208;
    1, psai(4, 1, 0, 1), -0.208;
    1, psai(4, 1, 1, 1), -0.208;
    2, psai(4, 0, 0, 2), -0.637;
    2, psai(4, 1, -1, 2), -0.291;
    2, psai(4, 1, 0, 2), -0.291;
    2, psai(4, 1, 1, 2), -0.291;
    };  
% 参与线性组合的原子轨道。一排一个轨道。
% 每排第一个数字表示轨道属于 R_List 中第几排坐标代表的原子，第二个元素为轨道，第三个为轨道能级，单位为Hartree（原子单位制下的能量单位）
% 原子轨道的具体形式请在 psai.m 中更改

filename = sprintf('GaAs_T，RMax=%g', RMax);  % 项目名称

if isempty(gcp('nocreate'))  % 预启动并行运算池
    parpool;
end

T = TMatrix(phi, R_List, R_L_List, potential); 
S = SMatrix(phi, R_List, R_L_List);  
H = HMatrix(phi, T, S); 
fprintf('Matrix generation complete! 喵。\n');

currentTime = datetime('now');
timeStr = datestr(currentTime, 'yyyy-mm-dd_HH-MM-SS');
name = sprintf('H, S_%s_%s', filename, timeStr);
save(name, 'H', 'S', 'filename', 'a', 'RMax')

% 运行main.m后，工作区中出现能量矩阵H和重叠矩阵S，之后可选择：
% 运行band_on_plane.m（负责计算）和band_on_plane_plot.m（负责绘图）以得到指定k平面上的能带图；
% 运行band_on_k_path以得到指定k路径上的能带图.