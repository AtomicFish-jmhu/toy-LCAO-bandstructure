function T = TMatrix(phi, R_List, R_L_List, potential)  % T矩阵
    syms x y z real;
    potential_total = 0;
    for k = 1:length(potential(:, 1))
        for j = 1:length(R_L_List(:, 1))
            potential_total = potential_total + subs(potential{k, 2}, [x, y, z], [ ...  % 总势场=所有原子的势场之和
                x-R_List(potential{k, 1}, 1)-R_L_List(j, 1), ...
                y-R_List(potential{k, 1}, 2)-R_L_List(j, 2), ...
                z-R_List(potential{k, 1}, 3)-R_L_List(j, 3), ...
            ]);
        end
    end

    N = length(phi(:, 1));
    T = sym(zeros(N));
    for m = 1:N
        for n = 1:N
            if n >= m
                fprintf('Computing T matrix element (%d, %d)...', m, n);
                T(m, n) = T_mn(phi(m, :), phi(n, :), R_List, R_L_List, potential_total, potential);
                if n == m
                    T(m, n) = real(T(m, n));
                end
            else
                T(m, n) = conj(T(n, m));  % 强制厄密化
            end
        end
    end
end

function T = T_mn(phi_m, phi_n, R_List, R_L_List, potential_total, potential)  % T矩阵元
    syms x y z R_x R_y R_z k_x k_y k_z real
    
    phi_n{2} = subs(phi_n{2}, [x, y, z], [x-R_x, y-R_y, z-R_z]);
    
    R_mn = [
        R_List(phi_n{1}, 1) - R_List(phi_m{1}, 1);
        R_List(phi_n{1}, 2) - R_List(phi_m{1}, 2);
        R_List(phi_n{1}, 3) - R_List(phi_m{1}, 3);
    ];

    R_L_num = length(R_L_List(:, 1));
    DeltaV_l = cell(R_L_num, 1);
    t_mn = zeros(R_L_num, 1);
    
    for i = 1:R_L_num
        V_l = subs(potential{phi_n{1}, 2}, [x, y, z], [ ...
            x-R_mn(1)-R_L_List(i, 1), ...
            y-R_mn(2)-R_L_List(i, 2), ...
            z-R_mn(3)-R_L_List(i, 3), ...
        ]);
        DeltaV_l{i} = potential_total - V_l;  % 减去R_L处原子的势
	DeltaV_l{i} = subs(DeltaV_l{i}, [x, y, z], [ ...  % 平移势场
		x + R_List(phi_m{1}, 1), ...
		y + R_List(phi_m{1}, 2), ...
		z + R_List(phi_m{1}, 3)  ...
		]);
    end
    
    integrand_List = cell(R_L_num, 1);  % 没法在并行循环体（parfor）里生成被积函数（涉及到符号运算），所以只能先生成好一个列表
    for i = 1:R_L_num
        integrand = conj(phi_m{2}) * DeltaV_l{i} * phi_n{2};
        integrand = matlabFunction(integrand, 'Vars', [x, y, z, R_x, R_y, R_z]); 
        integrand_List{i} = integrand;
    end
    
    parfor i = 1:R_L_num
        integrand = integrand_List{i};
        t_mn(i) = integral3( ...  % t_mn积分
            @(x, y, z) integrand( ...
                x, y, z, ...
                R_mn(1)+R_L_List(i,1), ...
                R_mn(2)+R_L_List(i,2), ...
                R_mn(3)+R_L_List(i,3)  ...
                ), ...
            -inf, inf, ...
            -inf, inf, ...
            -inf, inf, ...
            'AbsTol', 1e-3, ...
            'RelTol', 1e-2  ...
        );
    end 
    
    T = 0;
    for i = 1:R_L_num
        kdotsR = k_x*(R_mn(1) + R_L_List(i,1)) + ...  % 符号计算无法并行（在parfor中）进行，不然我就写在一起了
                 k_y*(R_mn(2) + R_L_List(i,2)) + ...
                 k_z*(R_mn(3) + R_L_List(i,3)); 
        T = T + exp(1i*kdotsR) * t_mn(i);
    end
    fprintf(' Done!\n')
end