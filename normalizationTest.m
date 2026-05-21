phi = {
    Y(0, 0);
    Y(1, -1);
    Y(1, 0);
    Y(1, 1);
    Y(2, -2);
    Y(2, -1);
    Y(2, 0);
    Y(2, 1);
    Y(2, 2);
    };

for i = 1:8
    vpa(norm(phi{i}))
end

function n = norm(Y_lm)
    syms x y z theta phi real;
    
    integrand = Y_lm^2;
    integrand = subs(integrand, [x, y, z], [ ...
        sin(theta) * cos(phi), ...
        sin(theta) * sin(phi), ...
        cos(theta) ...
        ]);
    integrand = integrand * sin(theta);
    
    n = int(int(integrand, theta, 0, pi), phi, 0, 2*pi);
    n = sqrt(n);
end

function angular = Y(l, m)  % 球对称势场下角向波函数为球函数，这里只写到d轨道，可以自行添加更多
    syms x y z real;
    r = sqrt(x^2 + y^2 + z^2);
    
    index = [l, m];
    switch l
        case 0
            N_lm = 1 / (2*sqrt(pi));
        case 1
            N_lm = (1/2) * sqrt(3/pi);
        case 2
            N_lm = (1/2) * sqrt(15/pi);
    end
    switch true
        case isequal(index, [0, 0])
            body = 1;
        case isequal(index, [1, -1])
            body = y / r;
        case isequal(index, [1, 0])
            body = z / r;
        case isequal(index, [1, 1])
            body = x / r;
        case isequal(index, [2, -2])
            body = x*y / r^2;
        case isequal(index, [2, -1])
            body = y*z / r^2;
        case isequal(index, [2, 0])
            body = (2*z^2 - x^2 - y^2) / (2 * sqrt(3) * r^2);
        case isequal(index, [2, 1])
            body = z*x / r^2;
        case isequal(index, [2, 2])
            body = (x^2 - y^2) / (2 * r^2);  
    end
    angular = N_lm * body;
end