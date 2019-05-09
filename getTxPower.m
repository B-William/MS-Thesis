%
% Function: getTxPower
% Inputs:
%   positions = (x, y) information for the users, given by
%   createNodeLayout()
%   SNR_dB = SNR of the weakest user in dB
%   Pnoise = noise power
%   pathLossCoeff = path loss coefficient (=2 for RF, =3 for certain UW)
% Outputs:
%   txPower = transmit power
%
function [ txPower ] = getTxPower( positions, SNR_dB, Pnoise, pathLossCoeff )

x0 = positions(1,1);
y0 = positions(1,2);
r = sqrt( (x0-positions(2:end,1)).^2 + (y0-positions(2:end,2)).^2 );
SNR_raw = 10^(SNR_dB/10);
rmax = max(r);
txPower = SNR_raw*Pnoise*rmax^pathLossCoeff;

end