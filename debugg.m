%% 我留下了这个脚本。用它来查看H，S在Gamma点的前2×2部分，用于调试

syms k_x k_y k_z;

H_dbg = subs(H, [k_x, k_y, k_z], [0, 0, 0]);
S_dbg = subs(S, [k_x, k_y, k_z], [0, 0, 0]);

vpa(H_dbg(1:2, 1:2))
vpa(S_dbg(1:2, 1:2))
