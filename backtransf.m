function [solution,time] = backtransf(x,y,tx,ep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BACKTRANSF transforms y from [0,1] to [ep,infty]
% FUNCTION CALL: [solution,time] = backtransf(x,y,tx,ep)
% INPUTS:  x    ... solution on [0,1]
%          y    ... solution on [0,1] which has to be transformed
%          tx   ... discrete time vector on [0,1]
%          ep   ... endpoint
% OUTPUTS: solution       
%          time
% AUTHOR: csimon
% DATE: 02/09
% COMMENT: none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tau2=sort(ep./tx(1:length(tx)));
if ep~= 0
     tau1=[];
     tau2=sort(ep./tx(1:length(tx)));
else    
    tau1=tx;
    tau2=sort(1./tx(1:length(tx)));
end 

tau=[tau1,tau2(2:end)]; %1 only once

y=fliplr(y);


solval=[x,y(2:end)];

solution=solval(1:end);
time=tau;


