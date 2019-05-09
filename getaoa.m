%%
% Function: getaoa
% Returns angle of arrival assuming tx node transmits to rx node
% Inputs:
%   tx = transmit node
%   rx = receive node
%   ula = linear array angles, given by createNodeLayout()
%   ang = angle matrix for the node layout, given by createNodeLayout()
% Outputs:
%   aoa = angle of arrival in degrees
%
function [aoa] = getaoa(tx,rx,ula,ang)
    aoa = mod(90 + ula(rx) - ang(rx,tx),180);
    if aoa > 90
        aoa = aoa - 180;
    elseif aoa < -90
        aoa = aoa + 180;
    end
end