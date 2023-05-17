clear
clc

n = 6;
k = 4;

x = [1:k-1, k-1];
s = 0;
i = 0;
while ~s
    i = i + 1;
    [s, x] = next_nchoosek(x, n, k);
end
assert(i == nchoosek(n, k), 'FAIL.');
fprintf('PASS.\n');
