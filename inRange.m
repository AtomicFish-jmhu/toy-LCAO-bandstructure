function R_L_List = inRange(RMax, T)  
    R_L_List = [];
    
    baseLength = sqrt(T(1, 1)^2 + T(1, 2)^2 + T(1, 3)^2);
    T_Normalized = T / baseLength;

    n = 10;  % 最大搜索范围
    for i = -n:n
        for j = -n:n
            for k = -n:n
                RL_Checking = [i, j, k];
                RL_Cartesian = RL_Checking * T_Normalized;
                RL_length = sqrt(RL_Cartesian(1)^2 + RL_Cartesian(2)^2 + RL_Cartesian(3)^2);
                if RL_length <= RMax + 0.1 
                    R_L_List = [R_L_List; RL_Checking];
                end
            end
        end
    end
    fprintf('共 %d 个原胞被包含在内。\n', length(R_L_List(:, 1)));
end