function [stopped, x] = next_nchoosek_bool(x, n, k)
% function [stopped, x] = next_nchoosek_bool(x, n, k)
% 
% Update the combination of choosing k values from integers 
% 1:n for one interation.
% 
% Input arguments:
%  x: boolean vector, whether each element in 1:n is 
%     choosen in current combination. It should be 
%     guarenteed that sum(x) == k.
%  n: Total number of objects in 1:n
%  k: Number of choosing objects from 1:n
% Output arguments:
%  stopped: 1 if the output "x" is the last combination, 
%           else 0
%  x: boolean vector, whether each element in 1:n is 
%     choosen in next combination. It should be guarenteed 
%     that sum(x) == k.

for i = 0:k-1
    if ~x(n-i)
        j = n-i;
        while ~x(j), j = j-1; end
        x(j:n) = false;
        x(j+1:j+1+i) = true;
        break
    end
end
stopped = all(x(1:n-k) == false(n-k,1));
return
