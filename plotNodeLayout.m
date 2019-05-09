%
% Function: plotNodeLayout
% Returns angle of arrival assuming tx node transmits to rx node
% Inputs:
%   positions = position information for each of the nodes
%   radius = simulation radius in meters
%   number = set to be 1 to number the nodes, set to 0 to not number
%
function plotNodeLayout( positions, radius, number )

N = size(positions,1);
xPositions = positions(:,1);
yPositions = positions(:,2);
labelshift = 10;
figure
scatter(xPositions(1),yPositions(1),50,'MarkerEdgeColor',[.9 .2 .2],'LineWidth',2)
hold on
scatter(xPositions(2:end),yPositions(2:end),50,'MarkerEdgeColor',[.2 .2 .7],'LineWidth',2)
if number
    text(xPositions+labelshift,yPositions+labelshift,cellstr(num2str([1:N]')),'FontSize',24)
end
grid on
xlabel('X Position (m)','FontSize',16)
ylabel('Y Position (m)','FontSize',16)
xlim([-radius radius])
ylim([-radius radius])
set(gca,'FontSize',24)
hold on

x = xPositions(1); y = yPositions(1);
r = radius;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit,yunit,'k--')
pbaspect([1 1 1])

end

