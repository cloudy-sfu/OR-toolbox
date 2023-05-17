clear
clc

m = 2; n = 5;

while 1
    A = randi(1000, [m, n]) ./ 100;
    b = randi(10000, [m, 1]) ./ 100;
    c = randn([1, n]) .* 50;
    
    [s1, x1] = revised_simplex(A, b, c, 'tol', 1e-13);
    [x2, ~, s2] = linprog(c', [],[], A, b, [0 0 0 0 0]);
    
    assert(s1 == s2, "Wrong fail code %d\n", s1);
    if s2 == 1
        err = norm(x1 - x2) / norm(x2);
        fprintf("Success success code %d\n", s2)
        fprintf("The error of solution %e\n", err)
    else
        fprintf("Correct fail code %d\n", s2);
    end
end