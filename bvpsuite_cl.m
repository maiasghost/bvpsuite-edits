function [outlog,optionsstruct]=bvpsuite_cl(bvpfile,optionsstruct);

% Idee: in optionsstruct sind alle zum Aufruf notwendigen Optionen
% gespeichert, ein Feld ist aber immer enthalten: optionsstruct.mode -
% damit wird bestimmt welche Subroutine aufgerufen wird, beispielsweise
% 'run', 'meshadaptation', etc.

% Base Mode if no other one is chosen: run;

if ~exist('optionsstruct')
    optionsstruct.mode='run';
end;

try
    
   if  strcmp(optionsstruct.mode,'run') | ...
       strcmp(optionsstruct.mode,'initialmesh') | ...
       strcmp(optionsstruct.mode,'meshadaptation') | ...
       strcmp(optionsstruct.mode,'equations') | ...
       strcmp(optionsstruct.mode,'errorestimate') 
   else       
       optionsstruct.mode = 'run';
   end
catch    
    optionsstruct.mode='run';
end

switch optionsstruct.mode
    case 'run'
        optionsstruct=randomarguments(bvpfile,optionsstruct);
        if strcmp(optionsstruct.start,'')
            if strcmp(optionsstruct.x1,'')
                optionsstruct.x1=feval(bvpfile,'x1');
            end
            ordnung=feval(bvpfile,'ordnung');
            [optionsstruct.x1,optionsstruct.start]=initialmesh(bvpfile,[optionsstruct.x1(1) optionsstruct.x1(end)],ones(size(ordnung,2),2),optionsstruct.x1,[],[],0);
        end        
        if strcmp(optionsstruct.plotrange,'')
            [outlog.coeff,outlog.x1,outlog.valx1,outlog.x1tau,outlog.valx1tau, ...
                outlog.polynomials,outlog.sol_infinite,outlog.tau_infinite, ...
                outlog.solx1_infinite,outlog.x1_infinite,outlog.lambda,outlog.eigenfunction]= ...
            run(bvpfile,optionsstruct.plot,optionsstruct.x1, ...
                optionsstruct.start,optionsstruct.bvpopt,optionsstruct.visible, ...
                optionsstruct.predictor);            
        else
            [outlog.coeff,outlog.x1,outlog.valx1,outlog.x1tau,outlog.valx1tau, ...
                outlog.polynomials,outlog.sol_infinite,outlog.tau_infinite, ...
                 outlog.solx1_infinite,outlog.x1_infinite,outlog.lambda,outlog.eigenfunction]= ...
            run(bvpfile,optionsstruct.plot,optionsstruct.x1, ...
                optionsstruct.start,optionsstruct.bvpopt,optionsstruct.visible, ...
                optionsstruct.predictor,optionsstruct.plotrange(1),optionsstruct.plotrange(2));
        end
    case 'initialmesh'
        optionsstruct=randomarguments(bvpfile,optionsstruct);
        [outlog.x1,outlog.coeff]= ...
        initialmesh(bvpfile,optionsstruct.x0, optionsstruct.valx0, ...
                    optionsstruct.x1, optionsstruct.par,optionsstruct.lambda,0);
%     case 'meshadaptation'
%         optionsstruct=randomarguments(bvpfile,optionsstruct);
%         [outlog.coeff,outlog.x1,outlog.valx1,outlog.x1tau,outlog.valx1tau, ...
%             outlog.polynomials,outlog.error1,outlog.tolerated]= ...
%         meshadaptation(optionsstruct.aTOL,optionsstruct.rTOL,optionsstruct.K, ...
%                         bvpfile,optionsstruct.plot,optionsstruct.x1, ...
%                         optionsstruct.start,optionsstruct.bvpopt, ...
%                         optionsstruct.visible,optionsstruct.maxiter, ...
%                         optionsstruct.finemesh,optionsstruct.predictor);
    case 'meshadaptation'
        optionsstruct=randomarguments(bvpfile,optionsstruct);
        if strcmp(optionsstruct.start,'')
            if strcmp(optionsstruct.x1,'')
                optionsstruct.x1=feval(bvpfile,'x1');
            end
            ordnung=feval(bvpfile,'ordnung');
            [optionsstruct.x1,optionsstruct.start]=initialmesh(bvpfile,[optionsstruct.x1(1) optionsstruct.x1(end)],ones(size(ordnung,2),2),optionsstruct.x1,[],[],0);
        end        
        if strcmp(optionsstruct.plotrange,'')
        [outlog.coeff,outlog.x1,outlog.valx1,outlog.x1tau,outlog.valx1tau, ...
            outlog.polynomials,outlog.error1,outlog.fine,outlog.sol_infinite,outlog.tau_infinite,outlog.solx1_infinite,outlog.x1_infinite,outlog.lambda,outlog.eigenfunction]= ...
         meshadaptation(optionsstruct.aTOL,optionsstruct.rTOL,optionsstruct.K, ...
                        bvpfile,optionsstruct.plot,optionsstruct.x1, ...
                        optionsstruct.start,optionsstruct.bvpopt, ...
                        optionsstruct.visible,optionsstruct.maxiter, ...
                        optionsstruct.finemesh,optionsstruct.predictor);
        else
        [outlog.coeff,outlog.x1,outlog.valx1,outlog.x1tau,outlog.valx1tau, ...
            outlog.polynomials,outlog.error1,outlog.fine,outlog.sol_infinite,outlog.tau_infinite,outlog.solx1_infinite,outlog.x1_infinite,outlog.lambda,outlog.eigenfunction]= ...
         meshadaptation(optionsstruct.aTOL,optionsstruct.rTOL,optionsstruct.K, ...
                        bvpfile,optionsstruct.plot,optionsstruct.x1, ...
                        optionsstruct.start,optionsstruct.bvpopt, ...
                        optionsstruct.visible,optionsstruct.maxiter, ...
                        optionsstruct.finemesh,optionsstruct.predictor,optionsstruct.plotrange(1),optionsstruct.plotrange(2));      
        end 
                    
    case 'equations'
        optionsstruct=randomarguments(bvpfile,optionsstruct)
        outlog.output = ...
        equations(optionsstruct.module,optionsstruct.coeff,bvpfile,optionsstruct.x1, ...
                  optionsstruct.outputpoints,optionsstruct.predictor,optionsstruct.psival, ...
                  optionsstruct.psi,optionsstruct.linear);
    case 'errorestimate'
        optionsstruct=randomarguments(bvpfile,optionsstruct);
        [outlog.x1tau,outlog.valerror,outlog.maxerrorvalx0,outlog.coeff_2, ...
            outlog.x1tau_2,outlog.val_2x1tau_2,outlog.x1_2,outlog.polynomials_2, ...
            outlog.val_2x1_2, outlog.valerror2,outlog.tau_infinite, ...
            outlog.error_infinite,outlog.error_infinite2,outlog.sol_infinite, ...
            outlog.lambda,outlog.eigenfunction]= ...
        errorestimate(optionsstruct.coeff,bvpfile,optionsstruct.plot, ...
                      optionsstruct.x1,optionsstruct.bvpopt,optionsstruct.visible, ...
                      optionsstruct.predictor);
%     case 'initialmesh2'
%         optionsstruct=randomarguments(bvpfile,optionsstruct);
%         [outlog.x1,outlog.start]= ...
%         initialmesh(bvpfile,optionsstruct.x0,optionsstruct.valx0, ...
%                     optionsstruct.x1,optionsstruct.par);
%     case 'pathfollowing'
%         optionsstruct=randomarguments(bvpfile,optionsstruct);
%         [outlog.coeff,outlog.x1,outlog.valx1,outlog.x1tau,outlog.valx1tau, ...
%             outlog.polynomials]= ...
%         pathfollowing(bvpfile,optionsstruct.plot,optionsstruct.x1, ...
%                       optionsstruct.start,optionsstruct.ersterschritt, ...
%                       optionsstruct.schrittanzahl,optionsstruct.bvpopt, ...
%                       optionsstruct.visible,optionsstruct.pfad, ...
%                       optionsstruct.startindex,optionsstruct.speichername, ...
%                       optionsstruct.pfaddatenhelp,optionsstruct.gitter, ...
%                       optionsstruct.finemesh);
    
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Function to assign random arguments %;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function optionsstruct=randomarguments(bvpfile,optionsstruct);

switch optionsstruct.mode
    case 'run'
        % plot - default 0;
        try
            if optionsstruct.plot == 1 | ...
               optionsstruct.plot == 0           
            else
                optionsstruct.plot=0;
            end
        catch
            15
            optionsstruct.plot=0;
        end
        
        % x1 - default '';
        try
            optionsstruct.x1;
        catch
            optionsstruct.x1='';
        end
        
        % start - default '';
        try 
            optionsstruct.start;
        catch
            optionsstruct.start='';
        end
        
        % bvpopt - default '';
        try
            optionsstruct.bvpopt;
        catch
            optionsstruct.bvpopt='';
        end
        
        %visible - default 0;
        try
            if optionsstruct.visible == 1 | ...
               optionsstruct.visible == 0
            else
               optionsstruct.visible=0;
            end
        catch
            optionsstruct.visible=0;
        end
        
        %predictor - default '';
        try
            optionsstruct.predictor;
        catch
            optionsstruct.predictor='';
        end
        
        try
            optionsstruct.plotrange;
        catch
            optionsstruct.plotrange='';
        end 
        
        
        
        
            
    case 'initialmesh'
        %x0 - default bvpfile-x1;
        try
            optionsstruct.x0;
        catch
            optionsstruct.x0=feval(bvpfile,'x1');
        end
        
        %valx0 - default 0;
        try
            optionsstruct.valx0;
        catch
            optionsstruct.valx0=zeros(length(feval(bvpfile,'ordnung')),length(optionsstruct.x0));
        end
            
       %x1 - default - bvpfile-x1;
       try
           optionsstruct.x1;
       catch
           optionsstruct.x1=feval(bvpfile,'x1');
       end
       
       %par - default - bvpfile-parameter;
       try
           optionsstruct.par;
       catch
           optionsstruct.par=feval(bvpfile,'parameter');
       end
       
       %lambda - default - [];
       try
           optionsstruct.lambda;
       catch
           optionsstruct.lambda=[];
       end
%     case 'meshadaptation'
%         % aTOL - default - 1e-6
%         try
%             optionsstruct.aTOL
%         catch
%             optionsstruct.aTOL=1e-6;
%         end
%         
%         % rTOL - default - 1e-6
%         try 
%             optionsstruct.rTOL
%         catch
%             optionsstruct.rTOL=1e-6;
%         end
%         
%         % K - default 100
%         try 
%             optionsstruct.K
%         catch
%             optionsstruct.K=100;
%         end
%         
%         %plot - default 0
%         try
%             if optionsstruct.plot == 0 | ...
%                optionsstruct.plot == 1
%             else
%                 optionsstruct.plot=0;
%             end
%         catch
%             optionsstruct.plot=0;
%         end
%         
%         %x1 - default linspace(a,b,50)
%         try
%             optionsstruct.x1
%         catch
%             a=feval(bvpfile,'x1');
%             optionsstruct.x1=linspace(a(1),a(end),50);
%         end
%         
%         %start - default ''
%         try
%             optionsstruct.start;
%         catch
%             optionsstruct.start='';
%         end;
%         
%         %bvpopt - default ''
%         try
%             optionsstruct.bvpopt;
%         catch
%             optionsstruct.bvpopt='';
%         end
%         
%         %visible - default 0
%         try
%             if optionsstruct.visible == 0 | ...
%                optionsstruct.visible == 1
%             else
%                 optionsstruct.visible=0;
%             end
%         catch
%             optionsstruct.visible=0;
%         end
%         
%         %maxiter - default 10
%         try
%             optionsstruct.maxiter;
%         catch
%             optionsstruct.maxiter=10;
%         end
%         
%         %finemesh - default 0
%         try
%             if optionsstruct.finemesh == 0 | ...
%                optionsstruct.finemesh == 1
%             else
%                 optionsstruct.finemesh=0;
%             end
%         catch
%             optionsstruct.finemesh=0;
%         end
%         
%         %predictor - default ''
%         try
%             optionsstruct.predictor;
%         catch
%             optionsstruct.predictor='';
%         end
%             
    case 'meshadaptation'
        % aTOL - default - 1e-6
        try
            optionsstruct.aTOL
        catch
            optionsstruct.aTOL=1e-6;
        end
        
        % rTOL - default - 1e-6
        try 
            optionsstruct.rTOL
        catch
            optionsstruct.rTOL=1e-6;
        end
        
        % K - default 100
        try 
            optionsstruct.K
        catch
            optionsstruct.K=100;
        end
        
        %plot - default 0
        try
            if optionsstruct.plot == 0 | ...
               optionsstruct.plot == 1
            else
                optionsstruct.plot=0;
            end
        catch
            optionsstruct.plot=0;
        end
        
        %x1 - default linspace(a,b,50)
        try
            optionsstruct.x1;
        catch
            a=feval(bvpfile,'x1');
            optionsstruct.x1=linspace(a(1),a(end),50);
        end
        
        %start - default ''
        try
            optionsstruct.start;
        catch
            optionsstruct.start='';
        end;
        
        %bvpopt - default ''
        try
            optionsstruct.bvpopt;
        catch
            optionsstruct.bvpopt='';
        end
        
        %visible - default 0
        try
            if optionsstruct.visible == 0 | ...
               optionsstruct.visible == 1
            else
                optionsstruct.visible=0;
            end
        catch
            optionsstruct.visible=0;
        end
        
        %maxiter - default 10
        try
            optionsstruct.maxiter;
        catch
            optionsstruct.maxiter=10;
        end
        
        %finemesh - default 0
        try
            if optionsstruct.finemesh == 0 | ...
               optionsstruct.finemesh == 1
            else
                optionsstruct.finemesh=0;
            end
        catch
            optionsstruct.finemesh=0;
        end
        
        %predictor - default ''
        try
            optionsstruct.predictor;
        catch
            optionsstruct.predictor='';
        end
    case 'equations'
        %was - default wert
        try
            optionsstruct.module;
        catch
            optionsstruct.module='value';
        end
        %if optionsstruct.module=='basispolynome0'
         %   optionsstruct.module='basispolynomials0';
        %end
        
        try 
            optionsstruct.coeff;
        catch
            optionsstruct.coeff='';
        end
        try 
            optionsstruct.x1;
        catch
            optionsstruct.x1=feval(bvpfile,'x1');
        end
        try
            optionsstruct.outputpoints;
        catch
            a=feval(bvpfile,'x1');
            optionsstruct.outputpoints=linspace(a(1),a(end),50);
        end
        try
            optionsstruct.predictor;
        catch
            optionsstruct.predictor='';
        end
        try 
            optionsstruct.psival;
        catch
            optionsstruct.psival='';
        end
        try
            optionsstruct.psi;
        catch
            optionsstruct.psi='';
        end
        try
            optionsstruct.linear;
        catch
            optionsstruct.linear='';
        end;           
    case 'errorestimate'
        %coeff - default '' - will not work - coeff is needed
        try
            optionsstruct.coeff;
        catch
            optionsstruct.coeff='';
        end
        
       %x1 - default ''
        try
            optionsstruct.x1;
        catch
            optionsstruct.x1='';
        end

        
        %plot - default 0
        try
            if optionsstruct.plot == 0 | ...
               optionsstruct.plot == 1
            else
                optionsstruct.plot=0;
            end
        catch
            optionsstruct.plot=0;
        end
        
        %bvpopt - default ''
        try
            optionsstruct.bvpopt;
        catch
            optionsstruct.bvpopt='';
        end
        
        %visible - default 0
        try
            if optionsstruct.visible == 0 | ...
               optionsstruct.visible == 1
            else
                optionsstruct.visible=0;
            end
        catch
            optionsstruct.visible=0;
        end

        %predictor - default ''
        try
            optionsstruct.predictor;
        catch
            optionsstruct.predictor='';
        end

       



end