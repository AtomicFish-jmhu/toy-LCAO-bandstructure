%% 我留下了这个脚本。用它来查看H，S，用于调试

syms k_x k_y k_z;

H_dbg = subs(H, [k_x, k_y, k_z], [0, 0, 0]);
S_dbg = subs(S, [k_x, k_y, k_z], [0, 0, 0]);

vpa(H_dbg)
vpa(S_dbg)
