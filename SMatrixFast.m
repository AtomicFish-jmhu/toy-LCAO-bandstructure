function S = SMatrixFast(phi, R_List, R_L_List)  % 重叠矩阵
    fprintf('SMatrix: Fast version deployed!\n');

    N = length(phi(:, 1));
    S = sym(zeros(N));
    for m = 1:N
        for n = 1:N
            if n >= m
                fprintf('Computing overlap matrix element (%d, %d)...', m, n);
                S(m, n) = S_mn(phi(m, :), phi(n, :), R_List, R_L_List);
                if n == m
                    S(m, n) = real(S(m, n));
                end
            else
                S(m, n) = conj(S(n, m));
            end
        end
    end
end

function S = S_mn(phi_m, phi_n, R_List, R_L_List)  % 重叠矩阵元
    syms x y z R_x R_y R_z k_x k_y k_z real
    
    phi_n{2} = subs(phi_n{2}, [x, y, z], [x-R_x, y-R_y, z-R_z]);
    
    integrand = conj(phi_m{2}) * phi_n{2};
    integrand = matlabFunction(integrand, 'Vars', [x, y, z, R_x, R_y, R_z]);
    
    R_mn = [
        R_List(phi_n{1}, 1) - R_List(phi_m{1}, 1);
        R_List(phi_n{1}, 2) - R_List(phi_m{1}, 2);
        R_List(phi_n{1}, 3) - R_List(phi_m{1}, 3);
    ];
    
    R_L_num = length(R_L_List(:, 1));
    s_mn = zeros(R_L_num, 1);
   
    parfor i = 1:R_L_num
        Rx = R_mn(1) + R_L_List(i,1);
        Ry = R_mn(2) + R_L_List(i,2);
        Rz = R_mn(3) + R_L_List(i,3);

        if Rx < 0.01 && Ry < 0.01 && Rz < 0.01
            d = 10;
        else
            d = 2 * sqrt(Rx^2 + Ry^2 + Rz^2) + 1e-15;  % 如果不加这个小量，采样器可能会算原子中心处的值，从而认为积分发散。但实际上被积函数具有伪奇异性
        end
        
        s_mn(i) = integral3(...  % 重叠积分
            @(x, y, z) integrand(x, y, z, Rx, Ry, Rz), ...
            Rx / 2 - d, Rx / 2 + d, ...
            Ry / 2 - d, Ry / 2 + d, ...
            Rz / 2 - d, Rz / 2 + d, ...
            'AbsTol', 1e-3, ... 
            'RelTol', 1e-2  ...
        );
    end
    
    S = 0;
    for i = 1:R_L_num
        kdotsR = k_x*(R_mn(1) + R_L_List(i,1)) + ...
                 k_y*(R_mn(2) + R_L_List(i,2)) + ...
                 k_z*(R_mn(3) + R_L_List(i,3)); 
        S = S + exp(1i*kdotsR) * s_mn(i);
    end
    fprintf(' Done!\n')
end