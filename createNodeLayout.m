%
% Function: createNodeLayout
% Inputs:
%   N = average # of nodes
%   radius = simulation layout radius
% Outputs:
%   positions = xy positions of each node
%   angles = relative angles between each node
%   ulaAng = random angular displacement of the linear array
%
function [ positions, angles, ulaAng ] = createNodeLayout( N, radius )

    res = 1e-10;
    theta = unifrnd(0,2*pi,N,1);
    r = radius*sqrt(unifrnd(0,1,N,1));
    xpos = r .* cos(theta);
    ypos = r .* sin(theta);
    positions = [0 0; xpos ypos];
        
    xpos = positions(:,1);
    ypos = positions(:,2);
    angles = zeros(N+1,N+1);    % row R, column C contains the angle of node C taking node R as the origin
    
    % calculates relative angles between each node
    for n = 1:N+1
        for nn = 1:N+1
            ydiff = ypos(nn)-ypos(n);
            xdiff = xpos(nn)-xpos(n);
            angles(n,nn) = atan2d(ydiff,xdiff);
            if angles(n,nn) < 0
                angles(n,nn) = angles(n,nn) + 360;
            end
        end
    end
    positions = round(positions/res)*res;
    ulaAng = randi([0 179],size(positions,1),1);
end

