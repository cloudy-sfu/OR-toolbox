clear
clc
load('revised_simplex_test_2.mat');

[s1, x1] = revised_simplex(A, b, c);
[x2, ~, s2] = linprog(c', [],[], A, b, [0 0 0 0 0]);

assert(s1 == s2);
if s2 == 1
    err = norm(x1 - x2) / norm(x2);
    disp(err);
else
    fprintf("Both returned fail code %d.\n", s2);
end
