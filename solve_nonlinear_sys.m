function [new_x,logstruct,itcount,fcount] = ...
    solve_nonlinear_sys(Fhandle,DFhandle,FDFhandle,x0,bvpfile,tau,bvpopt,parameters,psival,psi)
%   SOLVE_NONLINEAR_SYS  Solve nonlinear systems obtained by SBVPCOL and SBVPERR
%   
%   This routine is private and should not be accessed by other callers than SBVPCOL and SBVPERR   

% ********** store some fields of bvpopt as variables for quick reference
AbsTol      = bvpopt.AbsTol;
RelTol      = bvpopt.RelTol;
max_F_evals = bvpopt.ZfOpt.MaxFunEvals;
max_iter    = bvpopt.ZfOpt.MaxIter;
display     = bvpopt.ZfOpt.Display(1);
%display     = 'i';
Coutputfcn   = bvpopt.OutputFcn;
outputsel   = bvpopt.OutputSel;
trace       = bvpopt.OutputTrace;
%log         = bvpopt.Log;
log = true;
TRM         = bvpopt.TRM;

% ********** some initializations
fcount = 1;   % number of function evaluations
itcount= 1;   % number of iterations
new_x  = zeros(size(x0));
x = x0;
lambda = 1; 
lambdamin = bvpopt.Private.LambdaMin;
nPreviousTRMIterates = 0; 
logstruct = [];

if log
   logstruct.DOC = [];
   logstruct.tol_factor = [];
   logstruct.G0 = [];
   logstruct.G = [];
   logstruct.jac_update = [];
   logstruct.lambda = [];
   logstruct.nCorrections = [];
      
   logstruct.runtime.DFeval = [];
   logstruct.runtime.Feval = [];
   logstruct.runtime.lu = [];
   logstruct.runtime.resubstitute = [];
end

% Structure of the Algorithm:
% The Zerofinder consists of 3 different algorithms:
%    *) Fast Frozen Newton (cheap, but small domain of convergence (DOC) G1)
%    *) Predictor-Corrector Linesearch (expensive, but large DOC G2)
%    *) Trust Region Method (even more expensive, but even larger DOC G3)
%
% For the choice of the appropriate algorithm to use, it is necessary to determine the
% position of the current approximation w.r.t. the DOCs G1, G2, G3

% ********** Determine position of initial approximation

if display=='i'  % display information every iteration step
   fprintf('\n Determining zerofinder algorithm ... \n\n');
end


[TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
    determine_position(Fhandle,DFhandle,x,bvpfile,tau,bvpopt,parameters,psival,psi,lambdamin,fcount,logstruct,[],[]);   

if TolFactor < 1
   if display ~= 'o'  % off
      fprintf('\n\n Tolerances satisfied\n'); 
   end
   return
end

if DOC == 1
   x = new_x; % Accept new approximation
end

% ********** initialize Zerofinder Display
if display=='i'  % display information every iteration step
   switch DOC
       case 1 
           fprintf(' Fast Frozen Newton\n');
           fprintf(' STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');           
       case 2 
           fprintf(' Predictor-Corrector Line Search\n');
           fprintf(' STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
       case 3 
           fprintf(' Trust Region Method\n');
   end
end

% ********** Choice-of-Algorithm Loop
while 1       

if log logstruct.DOC = [logstruct.DOC DOC]; end

if itcount > max_iter
   %error(' Maximum number of iterations exceeded');
   return;
end

if fcount > max_F_evals
   %error(' Maximum number of function evaluations exceeded');
   return;
end

%Dasha Edit
% if fcount >= 5 || itcount >= 5
%     return;
% end

switch DOC

    
    % **************************************************************************   
    case 1 % ******************* Fast Frozen Newton ****************************
    % **************************************************************************   
      % variables that have to be available: a (=x), G0, G, F (=F(x)), simplified_delta_x
          
      % ********** Determine if Jacobian has to be updated    
      UpdateJac = G > bvpopt.Private.UpdateJacFactor * G0;
      
      if log 
         logstruct.lambda = [logstruct.lambda 1];
         logstruct.nCorrections = [logstruct.nCorrections 0];        
         logstruct.jac_update = [logstruct.jac_update UpdateJac];
      end
    
      % ********** Determine Jacobian
      if UpdateJac % Update Jacobian
                        cpt=cputime;
          DF = feval(DFhandle,bvpfile,x,tau,bvpopt,parameters,psival,psi);
                        if log 
                           logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                           logstruct.condestDF = condest(DF);
                        end

                        cpt=cputime;
          [L,U]=lu(DF);                  % LU-decomposition of DF
                        if log logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt]; end      

          lastwarn('');                  % initialize warning state            
            
                        cpt=cputime;
          delta_x = U\(L\(- F));         % Newton correction 
                        if log logstruct.runtime.resubstitute = ...
                           [logstruct.runtime.resubstitute cputime-cpt]; end      

          if length(lastwarn)            % Exit if DF is singular
             error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
          end    
            
          G0 = norm(delta_x);            % Norm of the Newton correction
      else % Iterate with frozen Jacobian
          delta_x = simplified_delta_x;
          G0 = G; % = norm(simplified_delta_x)
      end
      
      % ********** Perform Newton step
      new_x(:) = x(:) + delta_x;      % new approximation

                  cpt = cputime;           
      new_F = feval(Fhandle,bvpfile,new_x,tau,bvpopt,parameters,psival,psi);   % new residual
                  if log logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; end
      
      fcount = fcount +1;             % increase number of function evaluations

      % ********** Determine Position of new solution approximation
                  cpt = cputime;                
      simplified_delta_x = U\(L\(-new_F));            % simplified delta_x
                  if log logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt]; end      
      
      G = norm(simplified_delta_x);
      
      if log
         logstruct.G0 = [logstruct.G0 G0];
         logstruct.G =  [logstruct.G G];      
      end
      
      if G < (1-lambdamin/2) * G0 % Approximation improved
      
         % ***** Check if Tolerances are satisfied
         % ***** TolFactor: Factor by which the correction had to be
         % ***** reduced in order to satisfy the tolerances
         TolFactor = check_tolerances(new_x,simplified_delta_x,AbsTol,RelTol);
         
         if log logstruct.tol_factor = [logstruct.tol_factor TolFactor]; end
         
         if TolFactor < 1
             % Perform the last step towards the solution
             new_x(:) = new_x(:) + simplified_delta_x;
             
             if display=='i'  % display information every iteration step
                fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                  itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
             end

             if display ~= 'o'  % off
                fprintf('\n\n Tolerances satisfied\n'); 
             end

             return
         end

         DOC = 1;           % Keep iterating with the Fast Frozen Newton
         
         x = new_x;         % Accept new approximation
         F = new_F;         
         
         if display=='i'  % display information every iteration step
            fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                  itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
         end     
      elseif UpdateJac == 0 % Jacobian has not been updated in the previous step
         DOC = 1;           % -> Try FFN once again,
                            %    but update Jacobian this time
                            %    (We would have to evaluate the Jacobian
                            %    anyway in the Line Search Algorithm)

         if display=='i'  % display information every iteration step
            fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                  itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
         end
      else        
         DOC = 2;           % -> Switch to Predictor Corrector Line Search
         
         % ********** Prepare for next step with the respective method
         lambda = 1;
         
         % Use Jacobian computed here
         if log logstruct.jac_update  = [logstruct.jac_update  0]; end  
          
         if display=='i'
             fprintf(' %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                  itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);

             fprintf('\n Switching to Predictor-Corrector Line Search\n');
             fprintf(' STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
         end
      end   
                      
      itcount = itcount + 1;


   % **************************************************************************   
   case 2 % *************** Predictor-Corrector Line Search *******************
   % **************************************************************************   
   % Necessary variables at this point: a (=x), delta_x, simplified_delta_x, 
   % lambda (=lambda_pred), G0, G (=G(lambda_pred))   
   
   % Let G(lambda) = ||DF^(-1) * F(x + lambda*delta_x)|| 
   % At this point, x, delta_x, lambda, G(0), G'(0), G(lambda) are known
   % We check if lambda is accepted and correct lambda otherwise
      
   Accept_Lambda = G < (1-lambdamin/2)*G0;
   nCorrections = 0;

   while ~Accept_Lambda 
   % ********** Repeatedly correct lambda until it is accepted or a termination condition is met   
   
      if nCorrections == 0 % First correction, quadratical method
         % We know G(0), G'(0) and G=G(lambda) by previous computations.
         % To obtain a correction, we interpolate a quadratic polynomial through 
         % these values and define the minimum of this polynomial as the corrected lambda     
         % The respective polynomial reads G(x)=G(0) + x * G'(0) + x^2/lambda^2 *
         % (G(lambda) - G(0) - lambda * G'(0)). Zeroing the first derivative yields:      
   
         % store lambda (=lambda_pred) for possible further correction
         l2 = lambda;
         
         Gprime0 = -G0;
         lambda = - lambda * Gprime0 / (2* (G - G0 - lambda * Gprime0));
   
         % Allow lambda_cor between lambda_pred / 10 and 1
         lambda = max(l2/10 , min(lambda, 1));

         %***** Prepare for possible further correction (cubical, match notation with Num. Recipes)
         l1 = lambda;
                 
         nCorrections = 1;
      else % Subsequent corrections, cubical method
         coeff = 1/(l1-l2) * [1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2] * ...
                             [G1 - Gprime0*l1 - G0 ; G2 - Gprime0*l2 - G0];
            
         ca = coeff(1);
         cb = coeff(2);
                            
         lambda = (-cb + sqrt(cb^2-3*ca*Gprime0))/(3*ca);
             
         % Allow lambda_cor between lambda_cor_old / 10 and lambda_cor_old /2
         lambda = max(l1/10 , min(lambda, l1/2));
         
         % Prepare for possible further correction
         l2 = l1;
         l1 = lambda;
         
         nCorrections = nCorrections + 1;
      end    
          
      if lambda < lambdamin
         DOC = 3;  % switch to Trust Region Method 
         break;
      end
 
      new_x(:) = x(:) + lambda * delta_x;  % get new approximation
         
                  cpt = cputime;                   
      F = feval(Fhandle,bvpfile,new_x,tau,bvpopt,parameters,psival,psi);   % evaluate F at lambda_cor      
                  if log logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; end
      %Dasha Edit
      fcount = fcount + 1;
      fprintf('Function Count: %3i\n',fcount);  
                  cpt = cputime;                         
      simplified_delta_x = U\(L\-F);
                  if log logstruct.runtime.resubstitute = ...
                       [logstruct.runtime.resubstitute cputime-cpt]; end      
               
      G_old = G; % store old G for possible further correction
      G = norm(simplified_delta_x);
         
      Accept_Lambda = G < (1-lambdamin/2) * G0;
      
      % ***** Prepare for possible further correction
      G1 = G;
      G2 = G_old;
   end                          
            
   if Accept_Lambda % a lambda has been accepted
      x = new_x;    % Accept new solution
       
      if G < bvpopt.Private.SwitchToFFNFactor * G0    
         DOC = 1;
      else
         DOC = 2;
      end          
   end
   %Dasha Edit
   itcount = itcount + 1;
   fprintf('Iteration Count: %3i\n',itcount);
   if log 
      logstruct.lambda = [logstruct.lambda lambda]; 
      logstruct.nCorrections = [logstruct.nCorrections nCorrections];
      logstruct.G0 = [logstruct.G0 G0];
      logstruct.G =  [logstruct.G G];            
   end
   
   % At this point, DOC has been determined. Check Termination conditions if DOC ~= 3             
   if DOC ~= 3
      TolFactor = check_tolerances(x,simplified_delta_x,AbsTol,RelTol);
      
      if log logstruct.tol_factor = [logstruct.tol_factor TolFactor]; end
      
      if TolFactor < 1
         new_x(:) = x(:) + simplified_delta_x;
         
         if display=='i'  % display information every iteration step
            fprintf(' %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                    itcount,fcount,lambda,G/norm(x(:)),G0/G,TolFactor);
         end 

         if display ~= 'o'  % off
            fprintf('\n\n Tolerances satisfied\n'); 
         end
         
         return
      end
   elseif log
      logstruct.tol_factor = [logstruct.tol_factor -1];
   end

      

   
   if display=='i'  % display information every iteration step
      fprintf(' %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
              itcount,fcount,lambda,G/norm(x(:)),G0/G,TolFactor);
              
      if DOC == 1
          fprintf('\n Switching to Fast Frozen Newton\n');
          fprintf(' STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');           
      elseif DOC == 3
          fprintf('\n Switching to Trust Region Method\n');
      end
   end

   
   % If we keep up doing the Line Search, we have to predict lambda for the next step
   if DOC == 2   
          
      % ********** Keep old values of delta_x for prediction
      G0_old = G0;
          
      % ********** Determine delta_x
                        cpt = cputime;
      DF = feval(DFhandle,bvpfile,x,tau,bvpopt,parameters,psival,psi);
                        if log 
                           logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                           logstruct.jac_update = [logstruct.jac_update 1];
                           logstruct.condestDF = condest(DF);
                        end

                    cpt=cputime;
      [L,U]=lu(DF);                  % LU-decomposition of DF
                    if log logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt]; end      

      lastwarn('');                  % initialize warning state            
            
                    cpt=cputime;
      delta_x = U\(L\(- F));         % Newton correction 
                    if log logstruct.runtime.resubstitute = ...
                       [logstruct.runtime.resubstitute cputime-cpt]; end      

      if length(lastwarn)            % Exit if DF is singular
         error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
      end    

      G0 = norm(delta_x);
          
      % ********** Predict lambda [Deuflhardt74]
%dpr%      warning off MATLAB:divideByZero  % If Jacobians are the same, delta_x = simplified_delta_x
      mu = lambda * G0_old / norm(simplified_delta_x - delta_x);
%dpr%      warning on MATLAB:divideByZero 
         
      if mu > 0.7
         lambda = 1;
      else
         lambda = mu;
      end
          
      % ********** Make sure all necessary variables are available for correction step
      new_x(:) = x(:) + lambda * delta_x;
                    
                  cpt = cputime;                
      F = feval(Fhandle,bvpfile,new_x,tau,bvpopt,parameters,psival,psi);    
                  if log logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; end

      fcount = fcount + 1;
         
                  cpt = cputime;                
      simplified_delta_x = U\(L\-F);
                  if log logstruct.runtime.resubstitute = ...
                       [logstruct.runtime.resubstitute cputime-cpt]; end      

      G = norm(simplified_delta_x);                                       
   end

   
   % **************************************************************************   
   case 3 % ********************* Trust Region Method *************************
   % **************************************************************************   

   % Necessary variables at this point: a (=x), nPreviousTRMIterates
  %---------Achtung: Temporäre Änderung von mir!!!!
  if (~TRM)
      %Dasha Edit have uncommented this as not needed but returns function
      %output to whatever calls this
      err('trm');
      return;
  end
  %Ende Änderung
   % ********** Determine accuracy for Trust Region Method
   switch nPreviousTRMIterates
      case 0
         bvpopt.ZfOpt.TolX = sqrt(max(RelTol,AbsTol));
         bvpopt.ZfOpt.TolFun = 0;        % Assemble necessary options for FSOLVE
         bvpopt.ZfOpt.Jacobian = 'on';   % (Display, MaxIter, MaxFunEvals are user parameters)       
         bvpopt.ZfOpt.LargeScale = 'on';
      case 1
         bvpopt.ZfOpt.TolX = max(RelTol,AbsTol);
      case 2
         bvpopt.ZfOPt.TolX = 1000*eps;
      case 3
         error(' Requested tolerance is beyond the accuracy of the algorithm');
   end   
   % Dasha Edit
   bvpopt.ZfOpt.Display = 'iter';
   % *********** Limit the maximum number of iterations and FunEvals for TRM
   bvpopt.ZfOpt.MaxIter = max_iter - itcount;
   bvpopt.ZfOpt.MaxFunEvals = max_F_evals - fcount;                  
   
   % ********** Invoke Trust Region Method
   [x,dummy,F,ExitFlag,lsqLog,dummy2,DF] = ...
      lsqnonlin(FDFhandle,x,[],[],bvpopt.ZfOpt,bvpfile,tau,bvpopt,parameters,psival,psi);

   nPreviousTRMIterates = nPreviousTRMIterates +1;
   fcount = fcount + lsqLog.funcCount;
   icount = itcount + lsqLog.iterations;
  
   % ********** Intercept numerically singular Jacobians
   if length(findstr(lastwarn,'singular'))
       %dasha edit
      error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
      return;
   end    

   % ***** Determine whether anything has gone wrong
   %Dasha Edit lol
   if ExitFlag <= 0
     error(' Zerofinder did not converge. Increasing MaxIter and MaxFunEvals might help.');
     return;
    end
   
   
   % ********** Determine position of new approximation
   [TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
      determine_position(Fhandle,DFhandle,x,bvpfile,tau,bvpopt,parameters,psival,psi,lambdamin,fcount,logstruct,F,DF);   
   %dasha edit
   if TolFactor < 1 %1004.4 for length 2.5 100 points getting to 0.00625002 jacobian
      if display ~= 'o'  % off
         fprintf('\n\n Tolerances satisfied\n'); 
      end
      return
   end
 
   % ***** We only log those quantities if we continue iterating
   if log 
      logstruct.lambda = [logstruct.lambda -1];
      logstruct.nCorrections = [logstruct.nCorrections -1];      
      logstruct.jac_update = [logstruct.jac_update -1];
      logstruct.G0 = [logstruct.G0 G0];
      logstruct.G =  [logstruct.G G];
      logstruct.tol_factor = [logstruct.tol_factor TolFactor]; 
   end
  
 
   if DOC == 1
      x = new_x; % Accept new approximation
    elseif DOC == 2  % use Jacobian computed here
      if log logstruct.jac_update  = [logstruct.jac_update  0]; end %dpr%
   end
      
   
   if display=='i'  % display information every iteration step
      switch DOC
         case 1 
             fprintf(' Switching to Fast Frozen Newton\n\n');
             fprintf(' STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');           
         case 2 
             fprintf(' Switching to Predictor-Corrector Line Search\n\n');
             fprintf(' STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
         case 3 
             fprintf(' Continuing with Trust Region Method\n\n');
      end
   end
end
end



% *************************************************************************
% ********** DETERMINE POSITION OF INITIAL APPROXIMATION ******************
% *************************************************************************

function [TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
    determine_position(Fhandle,DFhandle,x,bvpfile,tau,bvpopt,parameters,psival,psi,lambdamin,fcount,logstruct,F,DF);   

log = bvpopt.Log;

% ********** Determine Newton correction
% ***** Evaluate Jabobian, if not provided (by fsolve)
if isempty(F)
             cpt = cputime;
   DF=feval(DFhandle,bvpfile,x,tau,bvpopt,parameters,psival,psi);
             if log 
                logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt]; 
                logstruct.condestDF = condest(DF);
             end      
             
             cpt = cputime;
   F = feval(Fhandle,bvpfile,x,tau,bvpopt,parameters,psival,psi);   % new residual
             if log logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; end                   

   fcount = fcount +1;            % increase number of function evaluations            
end

            cpt=cputime;
[L,U]=lu(DF);                  % LU-decomposition of DF
            if log logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt]; end      

lastwarn('');                  % initialize warning state            
            
            cpt=cputime;
delta_x = U\(L\(- F));         % Newton correction 
            if log logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt]; end      

if length(lastwarn)            % Exit if DF is singular
   error(' System matrix is close to singular. Try a refined mesh or another initial approximation.');
end    
            
G0 = norm(delta_x);            % Norm of the Newton correction

new_x = zeros(size(x));
new_x(:) = x(:) + delta_x;     % new approximation

% ***** In case we have very good initial approximation (e.g. from a previous mesh),
% ***** tolerances might be satisfied even now and we reduce computational
% ***** costs to the absolute minimum of one F- and one DF-evaluation
TolFactor = check_tolerances(new_x,delta_x,bvpopt.AbsTol,bvpopt.RelTol);

%Dasha Edit- end loop early
% if fcount == 20
%     TolFactor = 0.1;
%     return;
% end
if TolFactor < 1
   DOC = [];  % dummy outputs
   G = [];
   simplified_delta_x = [];
   
   if log
      logstruct.tol_factor = [logstruct.tol_factor TolFactor]; 
      logstruct.G0 = [logstruct.G0 G0];
   end
   
   return
end
   
% ********** Check Monotonicity Condition for lambda = 1 and lambda = lambda_min
            cpt = cputime;           
F = feval(Fhandle,bvpfile,new_x,tau,bvpopt,parameters,psival,psi);   % new residual
            if log logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; end
      
fcount = fcount +1;             % increase number of function evaluations

            cpt = cputime;                
simplified_delta_x = U\(L\(-F));            % simplified delta_x
            if log logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt]; end      
      
G = norm(simplified_delta_x);

if G <(1-lambdamin/2) * G0 % Approximation improved
   if G < bvpopt.Private.SwitchToFFNFactor * G0         % Approximation improved significantly
      DOC = 1;             % -> Fast Frozen Newton
      
      % See if tolerances are satisfied
      TolFactor = check_tolerances(new_x,simplified_delta_x,bvpopt.AbsTol,bvpopt.RelTol);
      if TolFactor < 1
         % Perform one last step towards the solution 
         new_x(:) = new_x(:) + simplified_delta_x;           

         if log
            logstruct.tol_factor = [logstruct.tol_factor TolFactor]; 
            logstruct.G0 = [logstruct.G0 G0];
            logstruct.G = [logstruct.G G];            
         end
         
         return
      end
   else % Approximation improved slightly
      DOC = 2;             % -> Predictor Corrector Line Search
      if log
         logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
      end
   end                
else % try another Newton step with smallest possible lambda
  
   new_x(:) = x(:) + lambdamin * delta_x;      % new approximation

            cpt = cputime;           
   Fmin = feval(Fhandle,bvpfile,new_x,tau,bvpopt,parameters,psival,psi);   % new residual
            if log logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt]; end
      
   fcount = fcount +1;             

            cpt = cputime;                
   simplified_delta_x_min = U\(L\(-Fmin));            
            if log logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt]; end      
      
   Gmin = norm(simplified_delta_x_min);
   
   if Gmin < (1-lambdamin/2) * G0 % Approximation improved
      DOC = 2; % Predictor Corrector Line Search
      if log 
         logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here      
      end
   else
      DOC = 3; % -> Trust Region Method
   end
end


%********************************************************************************
%******************************* CHECK TOLERANCES *******************************
%********************************************************************************

function TolFactor = check_tolerances(x,delta,AbsTol,RelTol)
% Calculates the factor by which delta has to be reduced to satisfy the
% tolerances. 
% out <= 1   => Tolerances satisfied
% out >  1   => Tolerances not satisfied
% ***** Determine Dimension of solution array x
[d p1 N] = size(x);

% ***** Extract solution components at mesh points (last meshpoint (b) not included
InfNorm_y = max(abs(reshape(x(:,1,:),d,N)));

% ***** Extract last correction at mesh points
delta = reshape(delta,size(x));
InfNorm_delta = max(abs(reshape(delta(:,1,:),d,N)));

% ***** Componentwise check if tolerances are satisfied
TolFactor = max(InfNorm_delta ./ (AbsTol + InfNorm_y * RelTol));
