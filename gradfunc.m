function ypr = gradfunc(x)
    n = length(x);
    ypr = zeros(n, 1);
    h = 1e-6;
    for i = 1:n
        x_forward = x;
        x_backward = x;
        x_forward(i) = x_forward(i) + h;
        x_backward(i) = x_backward(i) - h;
        ypr(i) = (func(x_forward) - func(x_backward)) / (2 * h);
    end
end
