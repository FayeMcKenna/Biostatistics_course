% this function receives a vector of p values - p - and a desired Q for FDR 
% it returns the corrected p_vlaue threshold that are believed by FDR (Q% of them, by
% mistake)
%
% [corrected_pval] = FDR(p, Q)
%

function [corrected_pval] = FDR(p, Q)

[sorted_pvals, sorted_idx] = sort(p);
n = length(p);
i = 1:n;
q = (Q * i/n)';
num_good = n;
for i = n:-1:1
    if q(i,1) < sorted_pvals(i)
        num_good = num_good - 1;
    else
        break  
    end
end 

% figure; 
% plot(sorted_pvals, 'b');
% hold on;
% title(ori);
% plot(q, 'r');
% X = 0:n/10:n;  % for the number of antigens in this experiment
% rectangle('Position', [0, Q, n, 1e-10]);

corrected_pval = max(sorted_pvals(1:num_good));