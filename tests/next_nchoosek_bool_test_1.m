clear
clc

n = 6;
k = 4;

x1(:,1) = 1:k;
x2 = zeros(n, 1);
x2(1:k, 1) = 1;
i = 1;
s1 = 0;
s2 = 0;
while ~s1
    i = i + 1;
    [s1, x1(:,end+1)] = next_nchoosek(x1(:,end), n, k);
    [s2, x2(:,end+1)] = next_nchoosek_bool(x2(:,end), n, k);
    assert(s1 == s2);
end
assert(i == nchoosek(n, k));

for j = 1:i
    assert(all( ...
        x1(:, j) == find(x2(:, j)) ...
    ));
end

disp("PASS.");
