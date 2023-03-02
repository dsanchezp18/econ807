% Model parameters
N = 10;          % number of grid points
delta = 1;       % depreciation rate
beta = 0.99;     % discount factor
alpha = 0.33;    % capital share
A = 1.8855;      % technology parameter

% Grid for capital accumulation
k_min = 0.1;     % minimum capital stock
k_max = 10;      % maximum capital stock
k_grid = linspace(k_min, k_max, N)';  % capital grid

% Utility function
u = @(c) log(c);  % log utility

% Compute closed-form solution
k_star = (alpha*beta*A)^(1/(1-alpha));
V_star = log(A*k_star^alpha - delta*k_star) / (1-beta);

% Initial guess for value function
V = zeros(N, 1);

% Value function iteration
tol = 1e-4;       % convergence tolerance
maxiter = 1000;   % maximum number of iterations
V_all = zeros(N, maxiter);  % value function at each iteration
max_V = zeros(maxiter, 1);  % maximum value function across iterations
for iter = 1:maxiter
    V_old = V;
    for i = 1:N
        k = k_grid(i);
        c = A*k^alpha + (1-delta)*k - k_grid;
        c(c<=0) = eps;
        u_vec = u(c);
        [V(i), ~] = max(u_vec + beta*V_old);
    end
    V_all(:, iter) = V;
    max_V(iter) = max(V);
    if norm(V - V_old) < tol
        break
    end
end

% Policy function
[~, k_idx] = max(u_vec + beta*V_old, [], 2);
k_policy = k_grid(k_idx);

% Plot value function iterations
figure
plot(k_grid, V_all(:, 1:iter))
title('Value Function Iterations')
xlabel('Capital Stock')
ylabel('Value')
legend('Iteration 1', 'Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', ...
    'Iteration 6', 'Iteration 7', 'Iteration 8', 'Iteration 9', 'Iteration 10', ...
    'Iteration 11', 'Iteration 12', 'Iteration 13', 'Iteration 14', 'Iteration 15')
saveas(gcf, 'value_function_iterations.png')

% Plot maximum value function
figure
plot(1:iter, max_V(1:iter))
title('Maximum Value Function')
xlabel('Iteration')
ylabel('Value')
saveas(gcf, 'maximum_value_function.png')

% Plot policy function and closed-form solution
figure
plot(k_grid, k_policy)
hold on
plot([k_min, k_max], [k_star, k_star], 'r--')
hold off
title('Policy Function')
xlabel('Capital Stock')
ylabel('Investment')
legend('Policy Function', 'Closed-Form Solution')
saveas(gcf, 'policy_function.png')

% Plot closed-form solution and final value function
figure
plot(k_grid, V_all(:, iter))
hold on
plot([k_min, k_max], [V_star, V_star], 'r--')
hold off
title('Value Function')
xlabel('Capital Stock')
ylabel('Value')
legend('Value Function', 'Closed-Form Solution')
saveas(gcf, 'value_function.png')

% --- With n = 100 ---- %

% Model parameters
N = 100;         % number of grid points
delta = 1;       % depreciation rate
beta = 0.99;     % discount factor
alpha = 0.33;    % capital share
A = 1.8855;      % technology parameter

% Grid for capital accumulation
k_min = 0.1;     % minimum capital stock
k_max = 10;      % maximum capital stock
k_grid = linspace(k_min, k_max, N)';  % capital grid

% Utility function
u = @(c) log(c);  % log utility

% Compute closed-form solution
k_star = (alpha*beta*A)^(1/(1-alpha));
V_star = log(A*k_star^alpha - delta*k_star) / (1-beta);

% Initial guess for value function
V = zeros(N, 1);

% Value function iteration
tol = 1e-4;       % convergence tolerance
maxiter = 1000;   % maximum number of iterations
V_all = zeros(N, maxiter);  % value function at each iteration
max_V = zeros(maxiter, 1);  % maximum value function across iterations
for iter = 1:maxiter
    V_old = V;
    for i = 1:N
        k = k_grid(i);
        c = A*k^alpha + (1-delta)*k - k_grid;
        c(c<=0) = eps;
        u_vec = u(c);
        [V(i), ~] = max(u_vec + beta*V_old);
    end
    V_all(:, iter) = V;
    max_V(iter) = max(V);
    if norm(V - V_old) < tol
        break
    end
end

% Policy function
[~, k_idx] = max(u_vec + beta*V_old, [], 2);
k_policy = k_grid(k_idx);

% Plot value function iterations
figure
plot(k_grid, V_all(:, 1:iter))
title('Value Function Iterations')
xlabel('Capital Stock')
ylabel('Value')
legend('Iteration 1', 'Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', ...
    'Iteration 6', 'Iteration 7', 'Iteration 8', 'Iteration 9', 'Iteration 10', ...
    'Iteration 11', 'Iteration 12', 'Iteration 13', 'Iteration 14', 'Iteration 15')
saveas(gcf, 'value_function_iterations_n100.png')

% Plot maximum value function
figure
plot(1:iter, max_V(1:iter))
title('Maximum Value Function')
xlabel('Iteration')
ylabel('Value')
saveas(gcf, 'maximum_value_function_n100.png')

% Plot policy function and closed-form solution
figure
plot(k_grid, k_policy)
hold on
plot([k_min, k_max], [k_star, k_star], 'r--')
hold off
title('Policy Function')
xlabel('Capital Stock')
ylabel('Investment')
legend('Policy Function', 'Closed-Form Solution')
saveas(gcf, 'policy_function_n100.png')

% Plot closed-form solution and final value function
figure
plot(k_grid, V_all(:, iter))
hold on
plot([k_min, k_max], [V_star, V_star], 'r--')
hold off
title('Value Function')
xlabel('Capital Stock')
ylabel('Value')
legend('Value Function', 'Closed-Form Solution')
saveas(gcf, 'value_function_n100.png')

% ---- With n = 1000 ---- %
% Model parameters
N = 1000;        % number of grid points
delta = 1;       % depreciation rate
beta = 0.99;     % discount factor
alpha = 0.33;    % capital share
A = 1.8855;      % technology parameter

% Grid for capital accumulation
k_min = 0.1;     % minimum capital stock
k_max = 10;      % maximum capital stock
k_grid = linspace(k_min, k_max, N)';  % capital grid

% Utility function
u = @(c) log(c);  % log utility

% Compute closed-form solution
k_star = (alpha*beta*A)^(1/(1-alpha));
V_star = log(A*k_star^alpha - delta*k_star) / (1-beta);

% Initial guess for value function
V = zeros(N, 1);

% Value function iteration
tol = 1e-4;       % convergence tolerance
maxiter = 1000;   % maximum number of iterations
V_all = zeros(N, maxiter);  % value function at each iteration
max_V = zeros(maxiter, 1);  % maximum value function across iterations
for iter = 1:maxiter
    V_old = V;
    for i = 1:N
        k = k_grid(i);
        c = A*k^alpha + (1-delta)*k - k_grid;
        c(c<=0) = eps;
        u_vec = u(c);
        [V(i), ~] = max(u_vec + beta*V_old);
    end
    V_all(:, iter) = V;
    max_V(iter) = max(V);
    if norm(V - V_old) < tol
        break
    end
end

% Policy function
[~, k_idx] = max(u_vec + beta*V_old, [], 2);
k_policy = k_grid(k_idx);

% Plot value function iterations
figure
plot(k_grid, V_all(:, 1:iter))
title('Value Function Iterations')
xlabel('Capital Stock')
ylabel('Value')
legend('Iteration 1', 'Iteration 2', 'Iteration 3', 'Iteration 4', 'Iteration 5', ...
    'Iteration 6', 'Iteration 7', 'Iteration 8', 'Iteration 9', 'Iteration 10', ...
    'Iteration 11', 'Iteration 12', 'Iteration 13', 'Iteration 14', 'Iteration 15')
saveas(gcf, 'value_function_iterations_n1000.png')

% Plot maximum value function
figure
plot(1:iter, max_V(1:iter))
title('Maximum Value Function')
xlabel('Iteration')
ylabel('Value')
saveas(gcf, 'maximum_value_function_n1000.png')

% Plot policy function and closed-form solution
figure
plot(k_grid, k_policy)
hold on
plot([k_min, k_max], [k_star, k_star], 'r--')
hold off
title('Policy Function')
xlabel('Capital Stock')
ylabel('Investment')
legend('Policy Function', 'Closed-Form Solution')
saveas(gcf, 'policy_function_n1000.png')

% Plot closed-form solution and final value function
figure
plot(k_grid, V_all(:, iter))
hold on
plot([k_min, k_max], [V_star, V_star], 'r--')
hold off
title('Value Function')
xlabel('Capital Stock')
ylabel('Value')
legend('Value Function', 'Closed-Form Solution')
saveas(gcf, 'value_function_n1000.png')

% --- New utility function, depreciation rate and also n = 100 --- %

% Model parameters
N = 100;         % number of grid points
delta = 0.025;   % depreciation rate
beta = 0.99;     % discount factor
alpha = 0.33;    % capital share
sigma = 2;       % sigma for the new utility function
A = 1.8855;      % technology parameter

% Grid for capital accumulation
k_min = 0.1;     % minimum capital stock
k_max = 10;      % maximum capital stock
k_grid = linspace(k_min, k_max, N)';  % capital grid

% Utility function
u = @(c) c.^(1-sigma)/(1-sigma);

% Value function iteration
v = zeros(N,1);  % initial guess for value function
tol = 1e-6;      % convergence tolerance
diff = Inf;      % initial difference
iter = 0;        % iteration counter
while diff > tol
    % Policy function
    c = (1-delta)*k_grid + A*k_grid.^alpha - k_grid(1);
    c(c<=0) = eps;  % ensure nonnegative consumption
    y = u(c) + beta*v;
    [v_new, k_idx] = max(y, [], 2);
    k = k_grid(k_idx);
    
    % Difference and update
    diff = max(abs(v_new - v));
    v = v_new;
    iter = iter + 1;
end

% Plot policy function
figure;
plot(k_grid, k, 'linewidth', 2);
title('Policy Function');
xlabel('Capital Stock');
ylabel('Next Period Capital Stock');

% Export plot as PNG
saveas(gcf, 'last_policy.png')
