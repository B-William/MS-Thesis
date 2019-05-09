%
% Function: chooseLNodes
% Inputs:
%   L = number of nodes to choose
%   range = range to choose over
%   txnode = node to exclude (tx-ing node)
% Outputs:
%   picks = L random nodes in range [min(range):max(range)] excluding txnode
%
function [picks] = chooseLNodes(L, range, txnode)

        list = min(range):max(range);
        n = length(list);

        % Choose L nodes that the txnode is also
        % transmitting to
        if n > L
            choices = randperm(n);
            picks = list(choices(1:L));
            if ~isempty(find(picks == txnode, 1))
               picks(picks == txnode) = choices(L+1);
            end
        else 
            picks = 1:n;
            picks(picks == txnode) = [];
        end
        picks = sort(picks);
end

