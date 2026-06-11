function [X, Y, Z] = band_structure3(H, S, baseLength, div, kz, R_List, R_L_List, phi)
    fprintf('Generating kpoints...\n');

    range = 6*pi/baseLength;
    
    kx = linspace(-range/2, range/2, div);
    ky = linspace(-range/2, range/2, div);
    [X, Y] = meshgrid(kx, ky);
    
    N = length(phi(:, 1));
    Z = zeros(div, div, N);
    
    R_L_num = length(R_L_List(:, 1));
    
    fprintf('Calculating bandstructure...');
    for i = 1:div
        parfor j = 1:div
            H_complete = zeros(N);
            S_complete = zeros(N);
            
            k_x = kx(i);
            k_y = ky(j);
            k_z = kz;
            
            for m = 1:N
                for n = 1:N
                    if n >= m
                        R_mn = [
                            R_List(phi{n, 1}, 1) - R_List(phi{m, 1}, 1);
                            R_List(phi{n, 1}, 2) - R_List(phi{m, 1}, 2);
                            R_List(phi{n, 1}, 3) - R_List(phi{m, 1}, 3);
                        ];
                        for k = 1:R_L_num
                            kdotsR = k_x*(R_mn(1) + R_L_List(k, 1)) + ...
                                     k_y*(R_mn(2) + R_L_List(k, 2)) + ...
                                     k_z*(R_mn(3) + R_L_List(k, 3)); 
                            H_complete(m, n) = H_complete(m, n) + exp(1i*kdotsR) * H(m, n, k);
                            S_complete(m, n) = S_complete(m, n) + exp(1i*kdotsR) * S(m, n, k);
                        end
                        if n == m
                            H_complete(m, n) = real(H_complete(m, n));
                            S_complete(m, n) = real(S_complete(m, n));
                        end
                    else
                        H_complete(m, n) = conj(H_complete(n, m)); 
                        S_complete(m, n) = conj(S_complete(n, m)); 
                    end
                end
            end
        
        solution = eig(H_complete, S_complete);
        Z(j, i, :) = real(solution);  % 注意 meshgrid() 返回的格式
        end
        if mod(i, floor(div * 0.1)) == 0
            fprintf('■');
        end
    end
    fprintf(' Complete!\n');
end