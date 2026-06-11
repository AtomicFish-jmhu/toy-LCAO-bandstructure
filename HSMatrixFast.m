function [H, S] = HSMatrixFast(R_List, R_L_List, potential, phi)
    fprintf('HSMatrix: Fast version deployed!\n')

    syms x y z R_x R_y R_z real

    potential_total = 0;
    for k = 1:length(potential(:, 1))
        for j = 1:length(R_L_List(:, 1))
            potential_total = potential_total + subs(potential{k, 2}, [x, y, z], [ ...  
                x-R_List(potential{k, 1}, 1)-R_L_List(j, 1), ...
                y-R_List(potential{k, 1}, 2)-R_L_List(j, 2), ...
                z-R_List(potential{k, 1}, 3)-R_L_List(j, 3), ...
            ]);
        end
    end  % 总势场=所有计入的原胞中所有原子的势场之和，于 R_List 中采用的绝对坐标下写出

    N = length(phi(:, 1));
    R_L_num = length(R_L_List(:, 1));
    H = zeros(N, N, R_L_num);  % 第三位存放 h_mn(R_l) （到计算能带时数值化 k 后再组装）
    S = zeros(N, N, R_L_num);  % 第三位存放 s_mn(R_l)

    for m = 1:N
        for n = 1:N
            if n >= m
                phi_n_shifted = subs(phi{n, 2}, [x, y, z], [x-R_x, y-R_y, z-R_z]); 
                
                R_mn = [
                    R_List(phi{n, 1}, 1) - R_List(phi{m, 1}, 1);
                    R_List(phi{n, 1}, 2) - R_List(phi{m, 1}, 2);
                    R_List(phi{n, 1}, 3) - R_List(phi{m, 1}, 3);
                ];
                
                fprintf('Computing Hamiltonian integrals of matrix element (%d, %d)...', m, n);
                    potential_local = subs(potential_total, [x, y, z], [ ... 
                        x + R_List(phi{m, 1}, 1), ...
                        y + R_List(phi{m, 1}, 2), ...
                        z + R_List(phi{m, 1}, 3)  ...
                        ]);  % 以 m 所在原子为原点

                    integrand = conj(phi{m, 2}) * Hamiltonian(phi_n_shifted, potential_local);
                    integrand = matlabFunction(integrand, 'Vars', [x, y, z, R_x, R_y, R_z]);

                    parfor i = 1:R_L_num
                        Rx = R_mn(1) + R_L_List(i, 1);
                        Ry = R_mn(2) + R_L_List(i, 2);
                        Rz = R_mn(3) + R_L_List(i, 3);

                        if abs(Rx) < 0.01 && abs(Ry) < 0.01 && abs(Rz) < 0.01
                            d = 10;
                        else
                            d = 2 * sqrt(Rx^2 + Ry^2 + Rz^2) + 1e-15; 
                        end

                        H(m, n, i) = integral3(...  
                            @(x, y, z) integrand(x, y, z, Rx, Ry, Rz), ...
                            Rx / 2 - d, Rx / 2 + d, ...
                            Ry / 2 - d, Ry / 2 + d, ...
                            Rz / 2 - d, Rz / 2 + d, ...
                            'AbsTol', 1e-3, ... 
                            'RelTol', 1e-2  ...
                        );
                    end
                    fprintf(' Done!\n')

                fprintf('Computing overlap integrals of matrix element (%d, %d)...', m, n);
                    integrand = conj(phi{m, 2}) * phi_n_shifted;
                    integrand = matlabFunction(integrand, 'Vars', [x, y, z, R_x, R_y, R_z]);

                    parfor i = 1:R_L_num
                        Rx = R_mn(1) + R_L_List(i, 1);
                        Ry = R_mn(2) + R_L_List(i, 2);
                        Rz = R_mn(3) + R_L_List(i, 3);

                        if abs(Rx) < 0.01 && abs(Ry) < 0.01 && abs(Rz) < 0.01
                            d = 10;
                        else
                            d = 2 * sqrt(Rx^2 + Ry^2 + Rz^2) + 1e-15;  
                        end

                        S(m, n, i) = integral3(...  
                            @(x, y, z) integrand(x, y, z, Rx, Ry, Rz), ...
                            Rx / 2 - d, Rx / 2 + d, ...
                            Ry / 2 - d, Ry / 2 + d, ...
                            Rz / 2 - d, Rz / 2 + d, ...
                            'AbsTol', 1e-3, ... 
                            'RelTol', 1e-2  ...
                        );
                    end
                    fprintf(' Done!\n')
            end
        end
    end
    fprintf('Matrix generation complete! 喵。\n');
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