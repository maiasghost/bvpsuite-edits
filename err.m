function err(type);

switch type
    %Errors from solve_nonlinear_sys.m
    case 'trm'
        error('The Newton solver did not converge! Possible reasons are: The starting guess for the solution values may be too imprecise; try to improve them! The initial mesh is too coarse; retry with smaller stepsize! The solution is not (locally) unique. Try other settings.');
    %Errors from equations.m
    case 'equations_err1'
        error('Point outside the defined interval!');
    case 'equations_err2'
        error('Value outside range');
    case 'equations_err3'
        error('Value outside range');
    case 'equations_err4'
        error('Number of points must be >=2 and <=9!');
    %Errors from bvpsuite.m
    case 'bvps_errdlg01'
        errordlg('Input a Matlab *.m file!','Error');
    case 'bvps_errdlg02'
        errordlg('You cannot use the name bvpsuite!','Error');
    case 'bvps_errdlg03'
        errordlg('The file does not exist. Fill in the form and klick save to make a new file or check the name of the file.','Error');
    case 'bvps_errdlg04'
        errordlg('The file last.log is missing.','Error');
    case 'bvps_errdlg05'
        errordlg('The file is incompatible with any version of bvpsuite!','Error');
    case 'bvps_errdlg06'
        errordlg('Input an m-file (Form: *.m)!','Error');
    case 'bvps_errdlg07'
        errordlg('You cannot overwrite bvpsuite itsself!','Error');
    case 'bvps_errdlg08'
        errordlg('The chosen file is not compatible with bvpsuite!','Error');
    case 'bvps_errdlg09'
        errordlg('Fill in the field Number / Partition rho_i!','Error');
    case 'bvps_errdlg10'
        errordlg('Gaussian points are only implemented with a maximum amount of 15. Choose "user" and determine them yourself.','Error');
    case 'bvps_errdlg11'
        errordlg('Avoid additional spaces in "Number / Partition rho_i"!','Error');
    case 'bvps_errdlg12'
        errordlg('Lobatto points are only implemented with a maximum amount of 15. Choose "user" and determine them yourself.','Error');
    case 'bvps_errdlg13'
        errordlg('Choose a value >= 2 for Lobatto points!','Error');
    case 'bvps_errdlg14'
        errordlg('Avoid additional spaces in "Number / Partition rho_i"!','Error');
    case 'bvps_errdlg15'
        errordlg('Bvpsuite can not cope with boundary conditions including derivatives at the right endpoint if the problem is posed on a semi-finite interval.','Error');
    case 'bvps_errdlg16'
        errordlg('Bvpsuite can not cope cope with boundary conditions including derivatives at the right endpoint if the problem is posed on a semi-finite interval.','Error');
    case 'bvps_errdlg17'
        errordlg('This file does not exist. Check the filename!','Error');
    case 'bvps_errdlg18'
        errordlg('The chosen file is not compatible with bvpsuite! Make sure that if your problem consists of one equation that the solution is denoted by z1 instead of z!.','Error');
    case 'bvps_errdlg19'
        errordlg('bvpsuite is not a bvpfile!','Error');
    case 'bvps_errdlg20'
        errordlg('For this option at least the field "lambda" is necessary!','Error');
    case 'bvps_errdlg21'
        errordlg('For this option the field "Initial mesh" is necessary!','Error');
    case 'bvps_errdlg22'
        errordlg('For this option the field "Initial values" is necessary!','Error');
    case 'bvps_errdlg23'
        errordlg('For this option the field "Alternative mesh" is necessary!','Error');
    case 'bvps_errdlg24'
        errordlg('The file "last.mat" does not exist, run a calculation to create it!');
    case 'bvps_errdlg25'
        errordlg('There are no inital values saved. Please fill in the field ''Initial values'' and save again.','Error');
    case 'bvps_errdlg26'
        errordlg('There is no inital mesh saved in this bvpfile','Error');
    case 'bvps_errdlg27'
        errordlg('Please fill in the field!','Error');
    case 'bvps_errdlg28'
        errordlg('The chosen file is not compatible with bvpsuite. The file was saved in an older version!','Error');
    case 'bvps_errdlg29'
        errordlg('The chosen file is not compatible with bvpsuite! Make sure that if your problem consists of one equation that the solution is denoted by z1 instead of z!.','Error');
    case 'bvps_errdlg30'
        errordlg('You need the symbolic math toolbox based on the Maple kernel to automatically create bvpfiles! This is extensively testet in Matlab Versions 7.0-7.2 (R14). Newer versions of Matlab (R2007b+) use the symbolic toolbox based on the MuPAD kernel which will unfortunately not work with bvpsuite. It is planned to implement the new kernel in future versions of bvpsuite, but at the moment you can only downgrade your Matlab version!','Error');
    case 'bvps_errdlg31'
        errordlg('Check your inputs! The error occurred when trying to convert your inputs in ''Equations'' to a Matlab inline-function, but it can also be an implication of inputs before (Orders, ...)!','Error');
    case 'bvps_errdlg32'
        errordlg('Check your inputs! The error occurred when trying to convert your inputs in ''Additiona-/Bondary Conditions'' to a Matlab inline-function, but it can also be an implication of inputs before (Orders, ...)!','Error');
    
    %error from meshadaptation.m
    case 'ma_maxiter'
        errordlg('The number of iterations allowed was succeeded - please check your settings and increase the maximum number of iterations');
    %error from initialmesh
    case 'initial1'
        errordlg('Insert at least two points in the interval [0,1].');
    case 'initial2'
        errordlg('Insert at least two points in the interval [1,L], where L is supposed to be large.');
        
    case 'version'
        errordlg('Using Matlab version 7.8. (R2009a) including the MuPAD based symbolic math toolbox an automatic transformation from a semi-infinite to a finite interval cannot be carried out. Use the Matlab versions 7.0-7.2 instead or transform the problem manually.');    
        
end;