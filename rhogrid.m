function grid = rhogrid(rho)
% Construct the staggered grid for the vetor function rho
% Written by GS, Innsbruck, 13 October 2005

% Staggered grid
M = length(rho);
offset = 1/M/2;
grid = linspace(offset,1-offset,M)';