%% 返回 H，S 在 Gamma 点的值

N = length(H(:, 1, 1));
R_L_num = length(H(1, 1, :));

k_x = 0;
k_y = 0;
k_z = 0;

H_complete = zeros(N);
S_complete = zeros(N);

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

% vpa(H_complete)
% vpa(S_complete)

