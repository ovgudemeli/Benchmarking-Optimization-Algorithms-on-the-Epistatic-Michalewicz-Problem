function y = func(x)

     n = length(x);
     m = 10;
     theta = pi / 6;
     y = 0;
    for i = 1:n
        if mod(i, 2) == 1 && i < n
            yi = x(i) * cos(theta) - x(i + 1) * sin(theta);
        elseif mod(i, 2) == 0 && i < n
            yi = x(i) * sin(theta) + x(i + 1) * cos(theta);
        else
            yi = x(i);
        end
        y = y - sin(yi) * (sin(i * yi^2 / pi))^(2 * m);
    end
end
  

