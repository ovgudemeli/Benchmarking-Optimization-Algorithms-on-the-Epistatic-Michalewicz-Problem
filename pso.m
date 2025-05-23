function pso(x)
    % Parameters
    n = 10;
    m = 10;
    theta = pi / 6;
    num_particles = 30;
    max_iter = 1000;
    
    % Constraints
    lb = zeros(1, n); % Lower bound
    ub = pi * ones(1, n); % Upper bound
    
    % PSO Parameters
    w = 0.729; % Inertia weight
    c1 = 1.49445; % Cognitive (particle)
    c2 = 1.49445; % Social (swarm)
    
    % Initial positions and velocities
    x = lb + (ub - lb) .* rand(num_particles, n);
    v = zeros(num_particles, n);
    pbest = x; % Best position of each particle
    pbest_val = arrayfun(@(i) objective(pbest(i, :), n, m, theta), 1:num_particles)';
    [gbest_val, gbest_idx] = min(pbest_val);
    gbest = pbest(gbest_idx, :);
    
    % Track the best solution value
    gbest_history = zeros(max_iter, 1);
    gbest_history(1) = gbest_val;
    
    % Optimization loop
    for iter = 2:max_iter
        % Update velocity and position
        for i = 1:num_particles
            r1 = rand(1, n);
            r2 = rand(1, n);
            v(i, :) = w * v(i, :) + c1 * r1 .* (pbest(i, :) - x(i, :)) + c2 * r2 .* (gbest - x(i, :));
            x(i, :) = x(i, :) + v(i, :);
            
            % Check bounds
            x(i, :) = max(min(x(i, :), ub), lb);
            
            % Evaluate
            val = objective(x(i, :), n, m, theta);
            if val < pbest_val(i)
                pbest(i, :) = x(i, :);
                pbest_val(i) = val;
            end
            if val < gbest_val
                gbest = x(i, :);
                gbest_val = val;
            end
        end
        
        % Save the best solution value
        gbest_history(iter) = gbest_val;
        
        % Iteration information
        if mod(iter, 100) == 0
            fprintf('Iteration %d: Best solution = %.4f\n', iter, gbest_val);
        end
    end
    
    % Results
    fprintf('Best solution = %.4f\n', gbest_val);
    fprintf('Best x = %s\n', mat2str(gbest));
    
    % Visualization
    figure;
    plot(1:max_iter, gbest_history, 'LineWidth', 2);
    title('Best Solution Over Iterations');
    xlabel('Iteration');
    ylabel('Best Solution Value');
    grid on;
end

function val = objective(x, n, m, theta)
    yi = zeros(1, n);
    for i = 1:n
        if mod(i, 2) == 1 && i ~= n
            yi(i) = x(i) * cos(theta) - x(i+1) * sin(theta);
        elseif mod(i, 2) == 0 && i ~= n
            yi(i) = x(i) * sin(theta) + x(i+1) * cos(theta);
        else
            yi(i) = x(i);
        end
    end
    val = -sum(arrayfun(@(i) sin(yi(i)) * (sin((i * yi(i)^2) / pi))^(2 * m), 1:n));
end
