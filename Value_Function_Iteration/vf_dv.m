% Housekeeping - Clear workspace
clear;
clc;


% Parameters
beta = 0.95;        % Discount factor
r = 0.03;           % Interest rate
sigma = 2;          % Coefficient of relative risk aversion
w = 1;              % Wage/income

% Discretization
a_min = 0;          % Minimum asset level
a_max = 100;        % Maximum asset level
num_a = 100;        % Number of asset grid points
a_grid = linspace(a_min, a_max, num_a);  % Asset grid

% Initialization of objects:
V = zeros(num_a, 1);        % Value function
V_new = zeros(num_a, 1);    % Updated value function
policy = zeros(num_a, 1);   % Policy function

% Utility function (CRRA)
u = @(c) (c^(1 - sigma) - 1) / (1 - sigma);

% Tolerance level and maximum iterations
tol = 1e-6;
max_iter = 1000;

% Value function iteration
for iter = 1:max_iter
    for i = 1:num_a
        a = a_grid(i);
        value = -Inf;
        for j = 1:num_a
            a_prime = a_grid(j);
            c = (1 + r) * a + w - a_prime;
            if c > 0
                value_new = u(c) + beta * V(j);
                if value_new > value
                    value = value_new;
                    V_new(i) = value;
                    policy(i) = c;
                end
            end
        end
    end
    
    % Check for convergence
    if max(abs(V_new - V)) < tol
        break;
    end
    
    V = V_new;
end

%%
% Display results:
disp('Value Function:')
disp(V)
disp('Policy Function (Consumption):')
disp(policy)

figure,plot(policy)



