%
% Function: plotArrays
% Inputs:
%   positions = position information for each user
%   ulaAng = angle of the linear array for each user in degrees
%   M = number of elements in the linear array
%   Fc = carrier frequency [Hz]
%
function plotArrays(positions, ulaAng, M, Fc)

radius = 100;
numUsers = size(positions,1);
d = (physconst('LightSpeed')/Fc)/2;
m = 0:M-1;
for user = 1:numUsers
    xy = positions(user,:).';
    pos = [d*m.' zeros(M,1)].'; pos(1,:) = pos(1,:) - (M-1)*d/2;
    ang = deg2rad(ulaAng(user));
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    rpos(:,:,user) = R*pos + xy;
end

figure
hold on
grid on
for user = 1:numUsers
    if user == 1
        col = [.9 .2 .2];
    else
        col = [.2 .2 .7];
    end
   scatter(rpos(1,:,user),rpos(2,:,user),80,col,'.')
   xPositions = rpos(1,1,user);
   yPositions = rpos(2,1,user);
   labelshift = 5;
   text(xPositions+labelshift,yPositions+labelshift,cellstr(num2str(user)),'FontSize',24)
end
xlabel('X Position (m)','FontSize',24)
ylabel('Y Position (m)','FontSize',24)
xlim([-radius radius])
ylim([-radius radius])
set(gca,'FontSize',24)
ticks = -100:50:100;
yticks(ticks)
xticks(ticks)
pbaspect([1 1 1])