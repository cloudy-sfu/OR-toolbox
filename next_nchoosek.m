function [stopped, x] = next_nchoosek(x, n, k)
% function [stopped, x] = next_nchoosek(x, n, k)
% 
% Update the combination of choosing k values from integers 
% 1:n for one interation.
% 
% Input arguments:
%  x: Current combination of k values from integers 1:n
%  n: Total number of objects in 1:n
%  k: Number of choosing objects from 1:n
% Output arguments:
%  stopped: 1 if the output "x" is the last combination, 
%           else 0
%  x: Next combination of k values from integers 1:n

for i = 0:k-1
    if x(k-i) < n-i
        x(k-i:k) = x(k-i) + (1:i+1);
        break
    end
end
stopped = x(1) == n - k + 1;
return
