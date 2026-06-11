function [X, Y, k_indices] = band_structure(H, S, kpoints, div, R_List, R_L_List, phi)
    fprintf('Generating dots sequence...');

    N = length(kpoints(:, 1)); 

    path_length = zeros(N-1, 1);
    for i = 1:(N-1)
        path_length(i) = sqrt((kpoints(i, 1) - kpoints(i+1, 1))^2 + ...
                              (kpoints(i, 2) - kpoints(i+1, 2))^2 + ...
                              (kpoints(i, 3) - kpoints(i+1, 3))^2);
    end

    points_num = zeros(N-1, 1);  
    total_points_num = 0;  
    for i = 1:(N-1)
        points_num(i) = ceil(path_length(i) * div);
        total_points_num = total_points_num + points_num(i);
    end
    total_points_num = total_points_num + 1; 

    kpath = zeros(total_points_num, 3); 
    trace = 0; 
    for i = 1:(N-1)
        for j = 1:3
            temp = linspace(kpoints(i, j), kpoints(i+1, j), points_num(i)+1).';
            kpath(trace + 1 : (trace + points_num(i)), j) = temp(1 : points_num(i));
        end
        trace = trace + points_num(i);
    end
    kpath(total_points_num, 1:3) = kpoints(N, 1:3); 

    k_indices = zeros(1, N);  
    k_indices(1) = 1;
    for i = 2:N
        k_indices(i) = k_indices(i-1) + points_num(i-1);
    end 

    fprintf('\nCalculating bands structure...');
    
    band_num = length(H(:, 1, 1));
    X = linspace(1, total_points_num, total_points_num);
    Y = zeros(total_points_num, band_num);

    R_L_num = length(R_L_List(:, 1));
    N = length(phi(:, 1));
    
    for i = 1:total_points_num
        H_complete = zeros(N);
        S_complete = zeros(N);
        
        k_x = kpath(i, 1);
        k_y = kpath(i, 2);
        k_z = kpath(i, 3);

        for m = 1:N
            for n = 1:N
                if n >= m
                    R_mn = [
                        R_List(phi{n, 1}, 1) - R_List(phi{m, 1}, 1);
                        R_List(phi{n, 1}, 2) - R_List(phi{m, 1}, 2);
                        R_List(phi{n, 1}, 3) - R_List(phi{m, 1}, 3);
                    ];
                    for j = 1:R_L_num
                        kdotsR = k_x*(R_mn(1) + R_L_List(j, 1)) + ...
                                 k_y*(R_mn(2) + R_L_List(j, 2)) + ...
                                 k_z*(R_mn(3) + R_L_List(j, 3)); 
                        H_complete(m, n) = H_complete(m, n) + exp(1i*kdotsR) * H(m, n, j);
                        S_complete(m, n) = S_complete(m, n) + exp(1i*kdotsR) * S(m, n, j);
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
        Y(i, :) = real(solution); 

        if mod(i, floor(total_points_num * 0.1)) == 0
            fprintf('■');
        end
    end
    fprintf(' Complete!\n');
end