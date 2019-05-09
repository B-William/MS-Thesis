%
% Function: makeCodes
% Inputs:
%   N = length of code
%   L = number of transmit beams
%   K = number of users
% Outputs:
%   code = matrix containing each users/beams code
%
function code = makeCodes(N, L, K)

if N > 1
    bc = 2*randi([0, 1], N, K)-1;
    bc_rep = repelem(bc, 1, L);
else
    bc = ones(N, K);
    bc_rep = repelem(bc, 1, L);
end

for i = 1:L:size(bc_rep, 2)
    rp = randperm(N);
    ind = rp(1:floor(N/2));
    bc_rep(ind, i) = -bc_rep(ind, i);
end
code = bc_rep;

end