function optionsstruct=options_cl(mode)

optionsstruct.mode=mode;

switch mode
    case 'run'
        % plot - default 0;
        optionsstruct.plot=0;
        optionsstruct.x1='';
        optionsstruct.start='';
        optionsstruct.bvpopt='';
        optionsstruct.visible=0;
        optionsstruct.plotrange='';
        
    case 'initialmesh'
        optionsstruct.x0='';
        optionsstruct.valx0='';
        optionsstruct.x1='';
        optionsstruct.par='';
        optionsstruct.lambda='';
      
    case 'meshadaptation'
        optionsstruct.aTOL=1e-6;
        optionsstruct.rTOL=1e-6;
        optionsstruct.K=100;
        optionsstruct.plot=0;
        optionsstruct.x1='';
        optionsstruct.start='';
        optionsstruct.bvpopt='';
        optionsstruct.visible=0;
        optionsstruct.maxiter=10;
        optionsstruct.finemesh=0;
        optionsstruct.plotrange='';
    case 'equations'
        optionsstruct.module='value';
        optionsstruct.coeff='';
        optionsstruct.x1='';
        optionsstruct.outputpoints='';
        optionsstruct.linear='';
    case 'errorestimate'
        optionsstruct.coeff='';
        optionsstruct.x1='';
        optionsstruct.plot=0;
        optionsstruct.bvpopt='';
        optionsstruct.visible=0;
    otherwise
        disp('No valid option chosen, please take a look at the manual.');
        optionsstruct='No valid option chosen, please take a look at the manual.';
     
end