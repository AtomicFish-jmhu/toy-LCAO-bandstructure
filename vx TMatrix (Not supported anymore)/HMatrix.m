function H = HMatrix(phi, T, S)  % 哈密顿矩阵
    N = length(phi(:, 1));
    H = sym(zeros(N));
    
    fprintf('Building Hamiltonian matrix...');
    for m = 1:N
        for n = 1:N
            H(m, n) = phi{n, 3} * S(m, n) + T(m, n);
        end
    end
    fprintf(' Done!\n')
end