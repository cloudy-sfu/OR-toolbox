clear
clc

load("biobjective_simplex_test_1.mat");

fprintf("=============== T2 ================\n");
[m, ~] = size(T2.A);
A_extended = [T2.A, eye(m)];
C_extended = [T2.C, zeros(2, m)];
[~, ~, lambda] = biobjective_simplex( ...
    A_extended, T2.b, C_extended, 'tol', 1e-15);
disp(lambda);

fprintf("=============== BLP 2017 ================\n");
[m, ~] = size(BLP2017.A);
A_extended = [BLP2017.A, eye(m)];
C_extended = [BLP2017.C, zeros(2, m)];
[~, ~, lambda] = biobjective_simplex( ...
    A_extended, BLP2017.b, C_extended, 'tol', 1e-15);
disp(lambda);

fprintf("=============== BLP 2018 ================\n");
[m, ~] = size(BLP2018.A);
A_extended = [BLP2018.A, eye(m)];
C_extended = [BLP2018.C, zeros(2, m)];
[~, ~, lambda] = biobjective_simplex( ...
    A_extended, BLP2018.b, C_extended, 'tol', 1e-15);
disp(lambda);
