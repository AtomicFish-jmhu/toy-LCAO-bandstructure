function S = SMatrix(phi, R_List, R_L_List)  % 重叠矩阵
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
        s_mn(i) = integral3(...  % 重叠积分
            @(x, y, z) integrand( ...
                x, y, z, ...
                R_mn(1) + R_L_List(i,1), ...
                R_mn(2) + R_L_List(i,2), ...
                R_mn(3) + R_L_List(i,3)  ...
            ), ...
            -inf, inf, ...
            -inf, inf, ...
            -inf, inf, ...
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