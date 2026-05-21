function H = HMatrix(phi, R_List, R_L_List, potential)  % 哈密顿矩阵
    syms x y z real;
    potential_total = 0;
    for k = 1:length(potential(:, 1))
        for j = 1:length(R_L_List(:, 1))
            potential_total = potential_total + subs(potential{k, 2}, [x, y, z], [ ...  % 总势场=所有原胞中所有原子的势场之和
                x-R_List(potential{k, 1}, 1)-R_L_List(j, 1), ...
                y-R_List(potential{k, 1}, 2)-R_L_List(j, 2), ...
                z-R_List(potential{k, 1}, 3)-R_L_List(j, 3), ...
            ]);
        end
    end

    N = length(phi(:, 1));
    H = sym(zeros(N));
    for m = 1:N
        for n = 1:N
            if n >= m
                fprintf('Computing Hamiltonian matrix element (%d, %d)...', m, n);
                H(m, n) = H_mn(phi(m, :), phi(n, :), R_List, R_L_List, potential_total);
                if n == m
                    H(m, n) = real(H(m, n));
                end
            else
                H(m, n) = conj(H(n, m));  % 强制厄密化
            end
        end
    end
end

function H = H_mn(phi_m, phi_n, R_List, R_L_List, potential_total)  % 哈密顿矩阵元/跃迁矩阵元
    syms x y z R_x R_y R_z k_x k_y k_z real
    
    phi_n{2} = subs(phi_n{2}, [x, y, z], [x-R_x, y-R_y, z-R_z]);  % 平移原子轨道
    
    potential_local = subs(potential_total, [x, y, z], [ ...  % 平移势场
        x + R_List(phi_m{1}, 1), ...
        y + R_List(phi_m{1}, 2), ...
        z + R_List(phi_m{1}, 3)  ...
        ]);
    
    integrand = conj(phi_m{2})*Hamiltonian(phi_n{2}, potential_local);
    integrand = matlabFunction(integrand, 'Vars', [x, y, z, R_x, R_y, R_z]);
    
    R_mn = [
        R_List(phi_n{1}, 1) - R_List(phi_m{1}, 1);
        R_List(phi_n{1}, 2) - R_List(phi_m{1}, 2);
        R_List(phi_n{1}, 3) - R_List(phi_m{1}, 3);
    ];
    
    R_L_num = length(R_L_List(:, 1));
    h_mn = zeros(R_L_num, 1);
    
    parfor i = 1:R_L_num
        Rx = R_mn(1) + R_L_List(i,1);
        Ry = R_mn(2) + R_L_List(i,2);
        Rz = R_mn(3) + R_L_List(i,3);
        
        h_mn(i) = integral3(...  % 跃迁积分
            @(x, y, z) integrand(x, y, z, Rx, Ry, Rz), ...
            -inf, inf, ...
            -inf, inf, ...
            -inf, inf, ...
            'AbsTol', 1e-3, ... 
            'RelTol', 1e-2  ...
        );
    end
    
    H = 0;
    for i = 1:R_L_num
        kdotsR = k_x*(R_mn(1) + R_L_List(i,1)) + ...  % 符号计算无法并行（在parfor中）进行，不然我就写在一起了
                 k_y*(R_mn(2) + R_L_List(i,2)) + ...
                 k_z*(R_mn(3) + R_L_List(i,3)); 
        H = H + exp(1i*kdotsR) * h_mn(i);
    end
    fprintf(' Done!\n')
end

function H = Hamiltonian(psi, potential) 
    syms x y z real
    hbar = 1.0;
    m = 1.0;
    
    laplacian = diff(psi, x, 2) + diff(psi, y, 2) + diff(psi, z, 2);
    kinetic = -hbar^2 / (2*m) * laplacian;
    
    potential = potential * psi;
    
    H = kinetic + potential;
end