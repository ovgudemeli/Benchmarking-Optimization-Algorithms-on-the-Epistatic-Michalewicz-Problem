function yprpr = hessianfunc(x)
    n = length(x);
    yprpr = zeros(n);
    h = 1e-6;
    for i = 1:n
        for j = 1:n
            x_ij1 = x; x_ij1(i) = x(i) + h; x_ij1(j) = x(j) + h;
            x_ij2 = x; x_ij2(i) = x(i) + h; x_ij2(j) = x(j) - h;
            x_ij3 = x; x_ij3(i) = x(i) - h; x_ij3(j) = x(j) + h;
            x_ij4 = x; x_ij4(i) = x(i) - h; x_ij4(j) = x(j) - h;
            yprpr(i,j) = (func(x_ij1) - func(x_ij2) - func(x_ij3) + func(x_ij4)) / (4 * h^2);
        end
    end
end