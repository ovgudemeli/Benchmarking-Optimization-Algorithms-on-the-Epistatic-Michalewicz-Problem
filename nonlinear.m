clc;
clear all;
close all;

realfmin= func([2.693,0.259,2.074,1.023,1.720]);
x0=pi*rand(1,5); % Başlangıç noktası
n=length(x0);
epsilon = 1e-5;
max_iter = 10;

% % Newton-Raphson
fprintf('\n--- Newton-Raphson ---\n');
x = x0';
path_newton = x';
tic;
for k = 1:max_iter
    g = gradfunc(x);
    H = hessianfunc(x);
    if rcond(H) < 1e-10 % Terslenebilirlik kontrolü
        warning('Hessian near-singular, terminating Newton-Raphson in iteration %d',k);
        break;
    end
    
    x_new = x - inv(H)* g;
    path_newton = [path_newton; x_new'];
    fprintf('Iter %d: f(x) = %.6f\n', k, func(x));
    if norm(g) < epsilon
        break;
    end
    x = x_new;
    
end
elapsed_time= toc;
 fvals_newton = arrayfun(@(i) func(path_newton(i,:)'), 1:size(path_newton,1));
    [min_val, idx] = min(fvals_newton);
    best_x = path_newton(idx,:);
    fprintf('Min value: %.4f at x = %s\n', min_val, mat2str(best_x, 4));
    fprintf("abs error= %.4f\n",abs(min_val-realfmin));
     fprintf('elapsed time is %f\n', elapsed_time);

% Hestenes-Stiefel
fprintf('\n--- Hestenes-Stiefel ---\n');
x = x0';
path_hs = x';
g = gradfunc(x);
d = -g;
tic;
for k = 1:max_iter
    alpha = 0.08; % Sabit adım büyüklüğü (isteğe göre ayarlanabilir)
    x_new = x + alpha * d;
    g_new = gradfunc(x_new);
    beta = (g_new' * (g_new - g)) / (d' * (g_new - g));
    d = -g_new + beta * d;
    path_hs = [path_hs; x_new'];
    fprintf('Iter %d: f(x) = %.6f\n', k, func(x));
    if norm(g_new) < epsilon
        break;
    end
    x = x_new;
    g = g_new;
  
end
elapsed_time=toc;
fvals_hs = arrayfun(@(i) func(path_hs(i,:)'), 1:size(path_hs,1));

[min_val, idx] = min(fvals_hs);
best_x = path_hs(idx,:);
fprintf('Min value: %.4f at x = %s\n', min_val, mat2str(best_x, 4));
fprintf("abs error= %.4f\n",abs(min_val-realfmin));
fprintf("elapsed time is %.6f\n", elapsed_time);

% % Polak-Ribiere
fprintf('\n--- Polak-Ribiere ---\n');
x = x0';
path_pr = x';
g = gradfunc(x);
d = -g;
tic;
for k = 1:max_iter
    alpha = 0.04;
    x_new = x + alpha * d;
    g_new = gradfunc(x_new);
    beta = (g_new' * (g_new - g)) / (g' * g);
    d = -g_new + beta * d;
    path_pr = [path_pr; x_new'];
    fprintf('Iter %d: f(x) = %.6f\n', k, func(x));
    if norm(g_new) < epsilon
        break;
    end
     x = x_new;
     g = g_new;
     
end
elapsed_time=toc;
    fvals_pr = arrayfun(@(i) func(path_pr(i,:)'), 1:size(path_pr,1));
    [min_val, idx] = min(fvals_pr);
    best_x = path_pr(idx,:);
    fprintf('Min value: %.4f at x = %s\n', min_val, mat2str(best_x, 4));
    fprintf("abs error= %.4f\n",abs(min_val-realfmin));
    fprintf("elapsed time is %.6f\n", elapsed_time);

%Fletcher-Reeves
fprintf('\n--- Fletcher-Reeves ---\n');
x = x0';
path_fr = x';
g = gradfunc(x);
d = -g;
tic;
for k = 1:max_iter
    alpha = 0.02;
    x_new = x + alpha * d;
    g_new = gradfunc(x_new);
    beta = (g_new' * g_new) / (g' * g);
    d = -g_new + beta * d;
    path_fr = [path_fr; x_new'];
    fprintf('Iter %d: f(x) = %.6f\n', k, func(x));
    if norm(g_new) < epsilon
        break;
    end
    x = x_new;
    g = g_new;
    elapsed_time=toc;
end
 fvals_fr = arrayfun(@(i) func(path_fr(i,:)'), 1:size(path_fr,1));
    [min_val, idx] = min(fvals_fr);
    best_x = path_fr(idx,:);
    fprintf('Min value: %.4f at x = %s\n', min_val, mat2str(best_x, 4));
    fprintf("abs error= %.4f\n",abs(min_val-realfmin));
    fprintf("elapsed time is %.6f\n", elapsed_time);


x1_vals = 0:0.05:pi;
x2_vals = 0:0.05:pi;
[X1, X2] = meshgrid(x1_vals, x2_vals);




% Sabit kalan x3-x10
x_fixed = 0;
Z = zeros(size(X1));

% Sadece x1 ve x2 değişiyor
for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        x_temp = x_fixed;
        x_temp(1) = X1(i, j);
        x_temp(2) = X2(i, j);
        Z(i, j) = func(x_temp');
    end
end
% grafiği 3B
figure;
surf(X1, X2, Z);
colorbar;
xlabel('x_1'); ylabel('x_2'); zlabel('f(x)');
title('Epistatic Michalewicz Function');

% 2B yüzey grafiği (kontur şeklinde)
figure;
contourf(X1, X2,Z, 20); % Dolgulu kontur
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Epistatic Michalewicz Function (x₁, x₂)');
grid on;


figure;
[C,h]=contourf(X1, X2,Z, 20); % Dolgulu kontur
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Newton Raphson (x₁, x₂)');
grid on;
hold on;
axis tight;
% Newton-Raphson yolu (kırmızı)
plot(path_newton(:,1), path_newton(:,2), 'r-o', 'LineWidth', 1, 'MarkerSize', 3);
uistack(h, 'bottom');

figure;
[C,h]=contourf(X1, X2,Z, 20); % Dolgulu kontur
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Hesteness Stiefel (x₁, x₂)');
grid on;
hold on;
axis tight;
axis([0 pi 0 pi]);
plot(path_hs(:,1), path_hs(:,2), 'b--s', 'LineWidth', 1, 'MarkerSize', 3);
uistack(h, 'bottom');

figure;
[C,h]=contourf(X1, X2,Z, 20); % Dolgulu kontur
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Polak Ribiere (x₁, x₂)');
grid on;
hold on;
axis tight;
axis([0 pi 0 pi]);
plot(path_pr(:,1), path_pr(:,2), 'g-.^', 'LineWidth', 1, 'MarkerSize', 3);
uistack(h, 'bottom');

figure;
[C,h]=contourf(X1, X2,Z, 20); % Dolgulu kontur
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Fletcher Reeves (x₁, x₂)');
grid on;
hold on;
axis tight;
axis([0 pi 0 pi]);
plot(path_fr(:,1), path_fr(:,2), 'm-d', 'LineWidth', 1, 'MarkerSize', 3);
uistack(h, 'bottom');
