%
% Function: getSMatrix
% Inputs:
%   spc = samples per chip
%   Ns = number of symbols in a packet
%   tau = each signals delay in samples
%   idletimes = each users idletimes in samples
%   code = matrix conttaining each signals code, code(:, 1) is signal one's
%   code
%   p = pulse shaping filter
%   users = list of active users/signals
% Outputs:
%   S = S matrix containing signature sequences at different delays
%   users = list of user indices corresponding to each column in S
%
function [S, users] = getSMatrix(spc, Ns, tau, idletimes, code, p, users)

    N = size(code,1);  
    ind = find(tau < 0);
    for i = 1:length(ind)
        if idletimes(ind(i)) <= spc*N*Ns
            tau(end+1) = idletimes(ind(i)); 
            code(:,end+1) = code(:,ind(i));
            users(end+1) = users(ind(i));
        end
    end
    
    K = length(tau);
    S = zeros(N*spc,K);
    S(1:spc:end,:) = code(:,:);
    S = repmat(S,[1 Ns]);
    if length(p) == spc
        S = filter(p,1,S);
    else
        S = filter(p,1,[S; zeros(floor(length(p)/2),size(S,2))]);
        S = S(floor(length(p)/2)+1:end,:);  
    end
    
    S = [zeros(abs(min(tau)),size(S,2)); S; zeros(max(tau)+(Ns-1)*spc*N,size(S,2))];
    for col = 0:size(S,2)-1
        symbolnumber = floor(col/K)+1;
        shift = tau(mod(col,K)+1) + (symbolnumber-1)*spc*N;
        S(:,col+1) = circshift(S(:,col+1),shift);
    end
    S = S/sqrt(spc*N); % normalize so that column inner product <s1, s1> = 1
end