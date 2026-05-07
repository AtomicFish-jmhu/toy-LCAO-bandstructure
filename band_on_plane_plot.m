%% 接收band_on_plane输出的X, Y, Z，绘制三维能带图。有三种画法可选

visual(X, Y, Z, 1)
% zlim([-100, 100]) 如果遇到发散困难

function visual(X, Y, Z, mode)  % 绘制指定k_x k_y平面上的能带图（在band_on_plane后使用）
    band_num = length(Z(1,1,:));
    if mode == 1  % 纵览所有点
        figure
        for k = 1:band_num
            plot3(X, Y, Z(:,:,k), '*', 'MarkerSize', 0.65)
            hold on;
        end
        hold off;
    end
    if mode == 2  % 查看各条能带（可能有部分点没有正确匹配，因此图像有瑕疵。在能带集中时尤甚）
        for k = 1:band_num
            figure(k)
            mesh(X, Y, Z(:,:,k))
        end
    end
    if mode == 3  % 指定截面图
        figure
        for k = 1:band_num
            plot(X(1,:), Z(:,1,k), '*', 'MarkerSize', 2, 'color', 'b')
            hold on;
        end
        hold off;
    end
end