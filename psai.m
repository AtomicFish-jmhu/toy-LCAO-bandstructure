function psi = psai(n, l, m, flag)  % matlab有个内置函数叫psi()，会重名，导致我没有用psi做函数名......
    psi = R(n, l, flag) * Y(l, m);
end

function radial = R(n, l, flag)  % 径向波函数
    syms x y z r real;
    
    if isequal([n, l, flag], [4, 0, 1])  % Ga, 4s
        zeta = 1.767;
    elseif isequal([n, l, flag], [4, 1, 1])  % Ga, 4p
        zeta = 1.555;
    elseif isequal([n, l, flag], [3, 2, 1])  % Ga, 3d
        zeta = 5.031;
    elseif isequal([n, l, flag], [4, 0, 2])  % As, 4s
        zeta = 2.236;
    elseif isequal([n, l, flag], [4, 1, 2])  % As, 4p
        zeta = 1.862;
    elseif isequal([n, l, flag], [3, 2, 2])  % As, 3d
        zeta = 5.793;
    else
        error('Not included!');
    end
    body = r^(n-1) * exp(-zeta*r);
    
    body_func = matlabFunction(body, 'Vars', r);
    N_nl = 1 / sqrt(integral(@(r) body_func(r).^2 .* r.^2, 0, Inf));  
    
    body = subs(body, r, sqrt(x^2 + y^2 + z^2));
    radial = N_nl * body;
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
            body = (2*z^2 - x^2 - y^2) / (2 * r^2);
        case isequal(index, [2, 1])
            body = z*x / r^2;
        case isequal(index, [2, 2])
            body = (x^2 - y^2) / (2 * r^2);  
    end
    angular = N_lm * body;
end