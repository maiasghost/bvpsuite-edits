function output = processV4(input)
% Process error data
% Written by GS, Innsbruck, 13 October 2005
% Rewritten Version V2, 16-18 October 2005, TU Wien.
% Version V3, 7-9 March, TU Wien
% Version V4, 21-22 Aug, TU Wien

% Rescale input function to have max = 1
scale = norm(input,inf); 
scaledinput = abs(input)/scale + 0.3;

% Limiter (lift small values away from zero)
scaledinput = limiterV3(scaledinput);

% Restore original scale
output = scale*scaledinput;

%============== Bias: see V3 ===================================== 