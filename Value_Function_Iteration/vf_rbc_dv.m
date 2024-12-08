% Value Function iteration example for a simple RBC model 
% Model has 2 state variables: capital and technology



% Housekeeping - Clear workspace
clear;
clc;

%% Parameters
beta = 0.95;        % Discount factor
alpha = 0.36;       % Capital share in production
delta = 0.08;       % Depreciation rate
rho = 0.9;          % Persistence of technology shock
sigma = 2;          % Coefficient of relative risk aversion
sigma_epsilon = 0.02;  % Std. dev. of technology shock


%% Discretization of objects
num_k = 50;         % Number of capital grid points
num_z = 6;          % Number of technology grid points
k_min = 0.1;        % Minimum capital level
k_max = 10;         % Maximum capital level
k_grid = linspace(k_min, k_max, num_k);  % Capital grid
[z_grid, z_prob] = rouwenhorst(num_z, 0, sigma_epsilon, rho); % Technology grid and transition matrix


%% Initialization
V = zeros(num_k, num_z);        % Value function
V_new = zeros(num_k, num_z);    % Updated value function
policy_k = zeros(num_k, num_z); % Policy function for capital

% Utility function (CRRA)
u = @(c) (c^(1 - sigma) - 1) / (1 - sigma);

% Tolerance level and maximum iterations
tol = 1e-6;
max_iter = 1500;


%% Value function iteration
for iter = 1:max_iter
    for i = 1:num_k
        k = k_grid(i);
        for j = 1:num_z
            z = z_grid(j);
            value = -Inf;
            for i_prime = 1:num_k
                k_prime = k_grid(i_prime);
                c = z * k^alpha + (1 - delta) * k - k_prime;
                if c > 0
                    expected_value = 0;
                    for j_prime = 1:num_z
                        expected_value = expected_value + z_prob(j, j_prime) * V(i_prime, j_prime);
                    end
                    value_new = u(c) + beta * expected_value;
                    if value_new > value
                        value = value_new;
                        V_new(i, j) = value;
                        policy_k(i, j) = k_prime;
                    end
                end
            end
        end
    end
    
    % Check for convergence
    if max(max(abs(V_new - V))) < tol
        break;
    end
    
    V = V_new;
end

%%
% Display results
disp('Value Function:')
disp(V)
disp('Policy Function (Capital):')
disp(policy_k)

figure,plot(policy_k)


