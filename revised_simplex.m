function [solved, x] = revised_simplex(A, b, c, max_iter, ...
    tol)
% function [solved, x] = revised_simplex(A, b, c, max_iter, ...
%     tol)
% Solving LP by revised simplex method.
% min   z = c*x
% s.t.  A*x = b,  b >= 0
%       x >= 0
% 
% Input arguments:
%  A: matrix, equality constraints
%  b: column vector, non-negative right-hand side
%  c: row vector, minimise objective function
%  max_iter: int, maximum iteration
%  tol: double, tolerance when calculating the reduced 
%       cost
% Returned values:
%  solved: status code (same as `linprog`)
%      1  linprog converged to a solution X.
%      0  Maximum number of iterations reached.
%     -2  No feasible point found.
%     -3  Problem is unbounded.
%  x: column vector, optimal feasible solution

arguments
    A (:,:) double
    b (:,1) double
    c (1,:) double
    max_iter (1,:) int32 = 500
    tol (1,:) double = 1e-15
end

% Validate arguments
[m, n] = size(A);
assert(length(b) == m, ['The length of constants b should ' ...
    'be equal to the number of rows of constraints A.']);
assert(length(c) == n, ['The length of costs c should be ' ...
    'equal to the number of columns of constraints A.']);
% Standard Computational form
assert(rank(A) == rank([A b]), ['Assumption should be ' ...
    'satisfied: rank(A) = rank([A b]).']);
assert(rank(A) == m, ['Assumption should be satisfied: ' ...
    'rank(A) = m.']);

x = zeros(n, 1);

% Basic feasible solution
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
reduced_cost = c - c(basis) / B * A;

% Update basis
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
    l_c = find(l_de > 0);  % candidate of λ
    % Problem is unbounded.
    if isempty(l_c), solved = -3; return; end
    [~, t] = min(l_nu(l_c) ./ l_de(l_c));
    t1 = find(basis);
    t = t1(l_c(t));
    
    basis(s) = true;
    basis(t) = false;
    B = A(:, basis);
    reduced_cost = c - c(basis) / B * A;
end

solved = 1;
x(basis) = B \ b;
end