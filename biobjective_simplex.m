function [solved, x, lambda] = biobjective_simplex( ...
    A, b, C, max_iter, tol)
% function [solved, x, lambda] = biobjective_simplex( ...
%     A, b, C, max_iter, tol)
%
% Solving bi-objective LP by revised simplex method.
% min   z = C*x
% s.t.  A*x = b,  b >= 0
%       x >= 0
% 
% Input arguments:
%  A: matrix, equality constraints
%  b: column vector, non-negative right-hand side
%  C: 2-row matrix, minimise objective function
%  max_iter: int, maximum iteration
%  tol: double, tolerance when calculating the reduced 
%       cost
% Returned values:
%  solved: status code (same as `linprog`)
%      1  linprog converged to a solution X.
%      0  Maximum number of iterations reached.
%     -2  No feasible point found.
%     -3  Problem is unbounded.
%  x: matrix, each column is an optimal feasible solution 
%     under the corresponding λ
%  lambda: row vector. It lists the thresholds when λ 
%          decreases from 1, the feasible basis of 
%          \bar{c} = λ*c1 + (1-λ)*c2 changes.

arguments
    A (:,:) double
    b (:,1) double
    C (2,:) double
    max_iter (1,:) int32 = 500
    tol (1,:) double = 1e-15
end


% Validate arguments
[m, n] = size(A);
assert(length(b) == m, ['The length of constants b should ' ...
    'be equal to the number of rows of constraints A.']);
assert(width(C) == n, ['The length of costs c should be ' ...
    'equal to the number of columns of constraints A.']);
% Standard Computational form
assert(rank(A) == rank([A b]), ['Assumption should be ' ...
    'satisfied: rank(A) = rank([A b]).']);
assert(rank(A) == m, ['Assumption should be satisfied: ' ...
    'rank(A) = m.']);

% Intialization
lambda = 1;
bar_c = C(1,:) * lambda(1) + C(2,:) * (1 - lambda(1));
x = zeros(n, 1);

% Phase 1: Basic feasible solution
basis = [true(m,1); false(n-m,1)];
B = A(:, basis);
x_b = B \ b;
stopped = 0;
while min(x_b) < 0  % feasible soltuion not found
    if stopped, solved = -2; return; end
    [stopped, basis] = next_nchoosek_bool(basis, n, m);
    B = A(:, basis);
    x_b = B \ b;
end
reduced_cost = bar_c - bar_c(basis) / B * A;

% Phase 2: Find optimal solution of LP(1)
I = eye(m);
iter = 0;
while min(reduced_cost) < - tol
    iter = iter + 1;

    % maximum number of iterations reached
    if iter >= max_iter, solved = 0; return; end

    % entering variable "s"
    [~, s] = min(reduced_cost(~basis));
    s1 = find(~basis);
    s = s1(s);
    
    % leaving variable "t"
    l_nu = I * (B \ b);  % numerator of λ
    l_de = I * (B \ A(:, s));  % denominator of λ
    l_c = find(l_de > 0);  % t_candidate
    % Problem is unbounded.
    if isempty(l_c), solved = -3; return; end
    [~, t] = min(l_nu(l_c) ./ l_de(l_c));
    t1 = find(basis);
    t = t1(l_c(t));
    
    % update basis
    basis(s) = true;
    basis(t) = false;
    B = A(:, basis);
    reduced_cost = bar_c - bar_c(basis) / B * A;
end
N = A(:, ~basis);

% Phase 3: Parametric biobjective simplex iterations
while lambda(end) > 0

    % entering variable "s"
    reduced_cost = C(:, ~basis) - C(:, basis) * ...
        (B \ N);
    bar_c1 = reduced_cost(1, :);
    bar_c2 = reduced_cost(2, :);
    % only negative element of \bar{c_N^2} are eligible
    t_c = find(bar_c2 < 0);
    if isempty(t_c), solved = 1; break; end
    % candidate of λ
    lambda_c = - bar_c2(t_c) ./ (bar_c1(t_c) - bar_c2(t_c));
    [lambda(end+1), s] = max(lambda_c);
    s1 = find(~basis);
    s = s1(t_c(s));
    
    % leaving variable "t"
    % It's an duplicate of finding leaving variables 
    % in Phase 2. I don't pack it as a new function because 
    % the input arguments are complex. It requires A, basis, 
    % and m, B derived from them. Re-calculating them in 
    % the function is low efficient.
    I = eye(m);
    r_nu = I * (B \ b);  % numerator of r
    r_de = I * (B \ A(:, s));  % denominator of r
    r_c = find(r_de > 0);  % r_candidate
    % Problem is unbounded.
    if isempty(r_c), solved = -3; return; end
    [~, t] = min(r_nu(r_c) ./ r_de(r_c));
    t1 = find(basis);
    t = t1(r_c(t));
    
    basis(s) = true;
    basis(t) = false;
    B = A(:, basis);
    N = A(:, ~basis);
    x(basis, end+1) = B \ b;
end

end