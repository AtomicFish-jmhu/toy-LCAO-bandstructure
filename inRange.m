function R_L_List = inRange(RMax, T)  % 搜索并添加截断半径RMax以内的格点，视原胞基矢里第一个的长度为1
    R_L_List = [];
    
    baseLength = sqrt(T(1, 1)^2 + T(1, 2)^2 + T(1, 3)^2);
    T_Normalized = T / baseLength;

    count = 0;
    n = 10;  % 这里只把搜索范围设定为10，事实上RMax超过5时已经很慢了（想马上出结果的话）......
    for i = -n:n
        for j = -n:n
            for k = -n:n
                RL_Checking = [i, j, k];
                RL_Cartesian = RL_Checking * T_Normalized;
                length = sqrt(RL_Cartesian(1)^2 + RL_Cartesian(2)^2 + RL_Cartesian(3)^2);
                if length <= RMax + 0.1
                    R_L_List = [R_L_List; RL_Checking];
                    count = count + 1;
                end
            end
        end
    end
%     fprintf('R_L_List = \n');
%     disp(R_L_List)
    fprintf('共 %d 个原胞被包含在内。\n', count);
end