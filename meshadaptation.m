function [coeff,x1,valx1,x1tau,valx1tau,polynomials,error1,sol_infinite,tau_infinite,lambda,eigenfunction,toleriert,ecol]=meshadaptation(aTOL,rTOL,K,bvpfile,plotsol,x1,start,bvpopt,plotres,maxiter,finemesh,predictor)


% meshadaptation ... calculate solution using adaptive meshes
%
% meshadaptation2(aTOL,rTOL,K,bvpfile,plotsol,x1,start,bvpopt,plotres,max
% iter) calculates the solution of a BVP, given in
% bvpfile using an adaptive mesh technique in order to solve the prescribed
% tolerance requirements given in aTOL (absolute tolerance) and rTOL
% (relative tolerance).
%
% K limits the maximal ratio between the shortest and the longest
% subinterval. plotsol (0,1) toggles the plot of the solution after the
% solution process. An initial mesh and an initial solution can be passed
% to the routine using x1 and start. If left empty, the mesh is determined
% from the bvpfile and the initial profile 0 is used. maxiter regulates the
% number of allowed adaptation steps.
%
% See also BVPSUITE, RUN
%


if ~exist('sol_infinite')
   
    sol_infinite=[];
end

if ~exist('tau_infinite')
   
    tau_infinite=[];
end

if ~exist('lambda')
   
    lambda=[];
end

if ~exist('eigenfunction')
   
  eigenfunction=[];
end
% prepare predictor
if ~exist('predictor') % predictor is set by the old bvpsuite routine - it serves no purpose here
    predictor=[];
end

% find initial mesh if not transmitted to meshadaptation2
if length(x1) < 1
    x1=feval(bvpfile,'x1');
end
left_end = x1(1);
right_end = x1(end);

% Calculate mesh density rho from the mesh
rho=(1./diff(x1)/length(x1))';
%Dasha Edit
%rhold = (1./diff(x1)/length(x1))';
%x1old = tau;
update_mode=1; % 1 .. update rho, 2 .. update N

M=length(x1)-1; % number of intervals
   
N = M-1;
Ncompare = 1e8; % initialised with a high value
standard=feval(bvpfile,'standard');
order_method = feval(bvpfile,'rho');
if standard
    order_method = order_method(2);
else
    order_method = length(order_method);
end

% Variables
minimprove = 0.1; % 0.1 means at least 10 % points less to try further density adjustment
safetyN = 0; % 0.1 means at least 10 % more points than actually suggested - not accurate anymore - use safety_sigma instead
safety_sigma = 0.9; % safety factor used in the formula of N-preestimation
LPFilter_k = order_method*2.75; % controls the gain of the adjustment of the mesh density - should be coppled to the method order
res_mode = 'tcolmid'; % triggers the choice of points for the evaluation of the residual

numequ=length(feval(bvpfile,'ordnung'));
if length(aTOL) <= numequ
    aTOL=ones(1,numequ)*aTOL(1);
end
if length(rTOL) <= numequ
    rTOL=ones(1,numequ)*rTOL(1);
end

parameter=feval(bvpfile,'parameter');

% DEFINITIONS ==================================
% Number of INTERVALS is M
% Number if internal GRID POINTS is N = M-1
% STEP SIZE on the auxiliary grid is dxi = 1/M
% rho is a vector with M components
% ==============================================

% Main loop: solution of problem + grid refinement
    
densityUpdates = 0;
% Dasha Edit- limit iterations
EarlyStopIteration = 0;
for j=1:maxiter
    EarlyStopIteration = EarlyStopIteration + 1;
    disp(['EarlyStopIteration' num2str(EarlyStopIteration)])
    % ***************************
    % Calculate next mesh from rho
    %*****************************
       
    % *************
    % SOLVE PROBLEM
    % *************
    
    [coeff,tau,y,tcol,ycol,polynomials] = run(bvpfile,0,x1,start,bvpopt,plotres,predictor);
    %Dasha Edit
%     rhold = (1./diff(x1)/length(x1))';
%     x1old = tau;   
    
    if j == 1
        if parameter>0
            par=coeff(length(coeff)-parameter+1:length(coeff));
        else
            par=[];
        end
    end

    if length(predictor)>0
        predictor=predictor(1:2,:);
    end
    % **************
    % ESTIMATE ERROR
    % **************
    [etcol,ecol,maxecol,ecoeff,tcol2,ycol2,tau2,polynomials2,y2,ecol2]=errorestimate(coeff,bvpfile,0,x1,bvpopt,plotres,predictor);    
  
    % Calculate error norms - global error ;
    %ycol_norm = max(abs(ycol),[],1)
    %err_norm = max(abs(ecol),[],1)
    ycol_norm = abs(ycol);
    err_norm = abs(ecol);
 
    %ycol_norm2 = max(abs(ycol2),[],1);
    %err_norm2 = max(abs(ecol2),[],1);
    ycol_norm2 = abs(ycol2);
    err_norm2 = abs(ecol2);

    for tolt=1:size(err_norm,1);
        [qTol(tolt),jmax(tolt)] = max(err_norm(tolt,:) ./ (aTOL(tolt) + ycol_norm(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
        [qTol2(tolt),jmax2(tolt)] = max(err_norm2(tolt,:) ./ (aTOL(tolt) + ycol_norm2(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
    end
    
    if plotres
        disp(['Tolerance Factor: ' num2str(qTol)])
    end
    if finemesh
        if plotres
            disp(['Tolerance Factor on the fine grid: ' num2str(qTol2)])
        end
    end
    %Dasha Edit, stops script early
%     if EarlyStopIteration == 2
%         qTol = 0.5;
%     end
    %%%%%%%%%%%%%%%%%%%%;
    % Update Procedure %;
    %%%%%%%%%%%%%%%%%%%%;
    
    if update_mode == 2; % update the number of mesh_points
        if qTol < 1
            if plotres
                disp('Tolerance satisfied')
            end
            x1=tau;
            valx1=y;
            x1tau=tcol;
            valx1tau=ycol;
            error1=0; % has no purpose ... not here, nor in meshadaptation!
            toleriert=0;
            break
        else
            if finemesh == 1
                if qTol2 < 1
                    if plotres
                        disp('Tolerance satisfied on the fine grid')
                    end
                    coeff=ecoeff;
                    x1=tau2;
                    valx1=y2;
                    x1tau=tcol2;
                    valx1tau=ycol2;
                    error1=0; % has no purpose ... not here, nor in meshadaptation!
                    toleriert=1;
                    break
                end
            end
            if (update_mode==2 && densityUpdates <2)
                update_mode = 1;
                densityUpdates = 1;
            else
                Nold = N;
                N = ceil((1+safetyN)*N*(max(max(err_norm,[],2)' ./ (safety_sigma*(aTOL))).^(1/(order_method+1))));
                %N = max(N,ceil(1.5^(1/order_method)*Nold)); % max increase
                % of number of points
            end
        end
    elseif update_mode == 1; % update the mesh density        
        if qTol < 1
            if plotres
                disp('Tolerance satisfied')
            end
            x1=tau;
            valx1=y;
            x1tau=tcol;
            valx1tau=ycol;
            error1=0; % has no purpose ... not here, nor in meshadaptation!
            toleriert=0;
            break
        else      
            if finemesh == 1
                if qTol2 < 1
                    if plotres
                        disp('Tolerance satisfied on the fine grid')
                    end
                    coeff=ecoeff;
                    x1=tau2;
                    valx1=y2;
                    x1tau=tcol2;
                    valx1tau=ycol2;
                    error1=0; % has no purpose ... not here, nor in meshadaptation!
                    toleriert=1;
                    break
                end
            end
            % Decide about the next step - Update N or rho ;
            Ncompare_old = Ncompare;
            Ncompare = ceil((1+safetyN)*N*(max(max(err_norm,[],2)' ./ (safety_sigma*(aTOL))).^(1/(order_method+1))));
            if update_mode == 1;
                if plotres
                    disp(['N suggested: ' num2str(Ncompare)])
                end
            end

            if Ncompare > (1-minimprove)*Ncompare_old % check if density is finally adjusted    
                % here we have to do the following - it can happen, that
                % the suggested number of points increases by a huge amount
                % so it would be a bad idea to update the density in this
                % last step - we have to use the older (better) density
                % here - this may happen in rare cases
                update_mode = 2;
                if Ncompare <= Ncompare_old        
                    N = Ncompare;
                    Ncompare = 1e8;
                else
                    N = Ncompare_old;
                    Ncompare = 1e8;
                    rho=rhold;
                    x1=x1old;
                end
            else
                densityUpdates = densityUpdates + 1;
            end
        end
    end
    
    % ******************
    % CALCULATE RESIDUAL
    % ******************
    if update_mode==1
        switch res_mode
        case 'mesh'
            resmesh = x1;
            residual=abs(equations('residual',coeff,bvpfile,tau,resmesh));
            % to get the residual for adjusting the mesh density we need a single
            % vector, so we take the 2-norm of the integral of it

            residualfinal=zeros(1,size(residual,2)-1);
            for i=1:size(residual,1)    
                residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(x1); % this means value on the lhs and rhs times the length
                residualfinal=residualfinal+residualadd.^2;
            end
            residualfinal=sqrt(residualfinal);
        case 'tcol'
            resmesh = tcol;
                       
            residual=abs(equations('residual',coeff,bvpfile,tau,resmesh));
            % to get the residual for adjusting the mesh density we need a single
            % vector, so we take the 2-norm of the integral of it

            residualfinal=zeros(1,length(x1)-1);
            for i=1:size(residual,1)    
                residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(tcol); 
                for k=1:length(x1)-1                    
                   residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)));
                end            
                residualfinal=residualfinal+residualadd2.^2;
            end
            residualfinal=sqrt(residualfinal);  
   
        case 'midpoints'
            resmesh = x1(1:end-1)+diff(x1)/2;
            residual=abs(equations('residual',coeff,bvpfile,tau,resmesh));
            residualfinal=zeros(1,size(residual,2));
            for i=1:size(residual,1)    
                residualadd=(residual(i,:)).*diff(x1); % this is the value at the midpoint times the length
                residualfinal=residualfinal+residualadd.^2;
            end
            residualfinal=sqrt(residualfinal);
        case 'tcolmid'
            resmesh = tcol(1:end-1)+diff(tcol)/2;
            residual=abs(equations('residual',coeff,bvpfile,tau,resmesh));
            % to get the residual for adjusting the mesh density we need a single
            % vector, so we take the 2-norm of the integral of it

            residualfinal=zeros(1,length(x1)-1);
            for i=1:size(residual,1)    
                residualadd=residual(i,:).*diff(tcol); 
                for k=1:length(x1)-1
                   residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)))/(x1(k+1)-x1(k));
                end            
                residualfinal=residualfinal+residualadd2.^2;
            end
            residualfinal=sqrt(residualfinal);
        end
    end
    
    Nnew = N;
    Mnew = Nnew+1;
    
    % Display stuff
    if plotres
        disp(['N used:      ' num2str(Nnew)])
        disp(['---------------'])
        disp(['Iteration ' num2str(j)])
    end
    if update_mode == 1
        if plotres
            disp(['Updating density'])
        end
        fprintf(1,'%s%d\n','Density Update N=',N);
    elseif update_mode == 2
        if plotres
            disp('Updating N')
        end
        fprintf(1,'%s%d\n','N Update N=',N);
    end
    rhold = rho;
    x1old = x1;
    if update_mode == 1;
        % Update density function
        rho = LPfilter(residualfinal',rho,LPFilter_k);   
    end;
 
    % check if the grid is within K-range
    maxrho = max(rho);
    rho=max(rho,maxrho/K);

    rho = bcTDFlogV4(rho);           % Smoothing
    rho = rho*sum(1./rho)/Mnew;      % Normalize
    
    x1 = (x1-x1(1))/(x1(end)-x1(1));
    I = cumtrapz(x1,([rho(1); rho]+[rho; rho(end)])'/2);
    x1 = pchip(I/I(end),x1,0:1/(Mnew):1);
    rho = 1./diff(x1)';
    x1 = left_end+(x1-0)*(right_end-left_end)/(1-0); %umrechnen von [0,1] auf [a,b]        
    [x1,start]=initialmesh2(bvpfile,tcol,ycol,x1,par);
    if length(predictor)>0
        predictor(3,1)=par(1);
    end
end 

% calculate eigenfunction, sol_infinite, lambda

 n=feval(bvpfile,'n');

    if feval(bvpfile,'Infinite')==1 &&  feval(bvpfile,'EVP')==1

        ep = feval(bvpfile,'Endpoint');

        if ep==0

            for i=1:n/2


                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);
                
                if i<n/2-1
                    eigenfunction(i,:)=sol_infinite(i,:);
                end 
                
            end

                
               % eigenfunction=sol_infinite(1:n/2-2,:);
                a=sol_infinite(n/2-1,1);
                lambda=a;
     
        else


            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

                if i<n-1
                    eigenfunction(i,:)=sol_infinite(i,:);
                end 
                
            end


            
                a=sol_infinite(n-1,1);
                lambda=a;
                

        end

    elseif feval(bvpfile,'Infinite')==1

        ep = feval(bvpfile,'Endpoint');

        if ep==0


            for i=1:n/2


                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);

            end

            
               
        else

            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end

        end




    elseif feval(bvpfile,'EVP')

     
            a=valx1tau(n-1,1);
            lambda=a;                       
            eigenfunction(1:n-2,:)=valx1tau(1:n-2,:);

    end
















% Draw the results if requested!
if plotsol
    figure;
    for i=1:length(ecol(:,1))
        subplot(length(ecol(:,1)),1,i);
        plot(etcol,real(ecol(i,:)),'color','black');
        axis tight
    end
    figure;
    plot(x1,ones(length(x1)),'.');
    axis tight
    figure;

    n=feval(bvpfile,'n');

    if feval(bvpfile,'Infinite')==1 &  feval(bvpfile,'EVP')==1



        ep = feval(bvpfile,'Endpoint');

        if ep==0

            for i=1:n/2


                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);

                %[valx1tau(i,:),x1tau]= backtransf(a(i,:),a(i+n/2,:),b);

            end


            for i=1:n/2-2
                subplot(n/2-2,1,i);
                plot(tau_infinite,-sol_infinite(i,:),'color','black');
                a=sol_infinite(n/2-1,1);
                title(['Eigenfunction for eigenvalue \lambda = ',num2str(a)])


            end

        else


            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end


            for i=1:n-2
                subplot(n-2,1,i);

                plot(tau_infinite,-sol_infinite(i,:),'color','black');
                a=sol_infinite(n-1,1);
                title(['Eigenfunction for eigenvalue \lambda = ',num2str(a)])
            end

        end





    elseif feval(bvpfile,'Infinite')==1

        ep = feval(bvpfile,'Endpoint');

        if ep==0


            for i=1:n/2


                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);

            end

            for i=1:n/2

                subplot(n/2,1,i);


                plot(tau_infinite,sol_infinite(i,:),'color','black');

            end
        else

            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end

            for i=1:n

                subplot(n,1,i);
                plot(tau_infinite,sol_infinite(i,:),'color','black');
                xlim([ep tau_infinite(end)])
            end






        end




    elseif feval(bvpfile,'EVP')

        for i=1:n-2
            subplot(n-2,1,i);
            plot(x1tau,valx1tau(i,:),'color','black');
            a=valx1tau(n-1,1);
            title(['Eigenfunction for eigenvalue \lambda = ',num2str(a)])


        end


    else

        for i=1:n
            subplot(n,1,i);
            plot(x1tau,valx1tau(i,:),'color','black');
            axis tight
        end

    end

end




if j==maxiter 
    %Dasha Edit
    err('ma_maxiter');
end



% **********************************************************
% **********************************************************
% Functions for processing the residual ********************
% **********************************************************
% **********************************************************


function tosignal = resampleV4(fromgrid,fromsignal,togrid)
% Oversampling that guarantees that tosignal remains positive
% Written by GS, TU Wien, 24 August 2006
% V4 is original version

signal = abs(fromsignal);                   % For robustness only
meansig = norm(signal,1)/length(signal);
mag = max(signal);
signal = signal/mag;                        % Normalize to maximum of 1
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-3
signal = log(signal);                       % Transform to logarithmic scale
hi = max(signal);
lo = min(signal);
%signal = pchip(fromgrid,signal,togrid); % Resample on new grid
signal = max(lo,min(hi,signal));         % Make sure amplitude doesn't grow
signal = exp(signal);                       % Now mapped to new grid and positive
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-2
signal = bcTDFlogV4(signal);                % Smoothing
signal = mag*signal/max(signal);            % Restore signal amplitude
newmean = norm(signal,1)/length(signal);    % Correct the signal's mean value
signal = signal + meansig - newmean;
signal = max(0,signal);                     % Protect against negative values
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-3
tosignal = signal;                          % Output

function newrho = LPfilter(err,rho,LPFilter_k)
% Generate step density profile update
% For use with 2pBVP solver for 1st order systems
% Euler-Lagrange optimal grid is generated using
% using deadbeat control law. Local error is
% equidistributed over the interval.

% err is already on staggered grid in this version
k = LPFilter_k;       % k = 3 for DBC, k = 4 conv filter (root 1/3)
M = length(rho);

% Process input error 
scalederr = M*abs(err).^(1/k);
% scalederr = processV4(bcTDFlogV4(scalederr)); % SWITCHED OFF

ratio = scalederr;
% ratio is suggested update of rho; compute next rho
rhonew = rho.*ratio;
%rhonew = bcTDFlogV4(rhonew);                   % SWITCHED OFF
newrho = rhonew;%*sum(1./rhonew)/M;


function out = bcTDFlogV4(in)
% Boundary corrected Topelitz digital filter, multiplicative (logarithmic)
% Applies 2nd order convolution LP filter repeatedly
% to input signal until output signal is sufficiently smooth

smooth = 50;     % Define smoothness requirement
N = length(in);
old = in;
signal = old;    % Mem alloc
snratio = 0;

while snratio < smooth
    % Apply filter until S/N ratio at this stage is at least = smooth 
    for i=2:N-1
        signal(i) = (old(i-1)*old(i)^2*old(i+1))^(1/4);
    end
    % Boundary corrections
    signal(1) = (old(1)^3*old(2)^3/old(3)^2/old(4)*old(5))^(1/4);
    % signal(1) = (old(1)^2*old(2)^3/old(3))^(1/4);
    signal(N) = (old(N)^3*old(N-1)^3/old(N-2)^2/old(N-3)*old(N-4))^(1/4);
    % signal(N) = (old(N)^2*old(N-1)^3/old(N-2))^(1/4);
    
% Compute noise    
    noise = old - signal;
    s = norm(signal);
    n = norm(noise);
    snratio = s/(n + 1e-3*s);
    old = signal;
end
out = signal;


function grid = rhogrid(rho)
% Construct the staggered grid for the vector function rho

% Staggered grid
M = length(rho);
offset = 1/M/2;
grid = linspace(offset,1-offset,M)';



% This one is used only here - later initialmesh should serve all purposes
function [x1,start]=initialmesh2(bvpfile,stellen,werte,x1,par)
if (nargin<6) koeff=[]; if (nargin<5) par=[]; end;end;

N=length(x1)-1;
n=feval(bvpfile,'n');
m=feval(bvpfile,'m');
ordnung=feval(bvpfile,'ordnung');
parameter=feval(bvpfile,'parameter');
dummystart=ones(N*(n*m+sum(ordnung))+parameter,1);
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];
if length(par)==0
    par=ones(parameter);
end
if length(koeff)>0 && parameter>0
    par=koeff(end-parameter+1:end);
end
if max(ordnung)>0
    psi=equations('basispolynome',dummystart,bvpfile,x1);
end
if min(ordnung)==0
    psi0=equations('basispolynome0',dummystart,bvpfile,x1);
end
for i=0:N-1
    h(i+1)=x1(i+2)-x1(i+1);
end


for ord=1:max(ordnung)
    proj_intervallstellen=0:1/(m+ord-1):1;
    for i=1:m
        psival(ord,i,1:m+ord)=polyval(psi(i,:,ord),proj_intervallstellen);
    end
end
jj=0;
for ii=2:length(stellen)
    if stellen(ii)~=stellen(ii-1)
        jj=jj+1;
        stellen2(jj)=stellen(ii-1);
        werte2(1:n,jj)=werte(1:n,jj);
    end
end
stellen2(jj+1)=stellen(end);
werte2(1:n,jj+1)=werte(1:n,end);
for p=1:n
    ordp=ordnung(p);
    interpolationsstellen=zeros(1,(N-1)*(ordnung(p)+m-1)+1);
    for i=0:N-1
        interpolationsstellen(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
    end
    interpolationsstellenwerte=interp1(stellen2,werte2(p,:),interpolationsstellen,'spline');
    for i=0:N-1
        if (m==1 && ordnung(p)==0)
            intervallstellen(1)=x1((i+1));
        else
            intervallstellen=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
        end
      % if length(koeff)==0
            %intervallwerte=interp1(stellen,werte(p,:),intervallstellen,'spline');
            intervallwerte=interpolationsstellenwerte(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1);
      %  else
      %      intervallwerte=equations('wert',koeff,bvpfile,stellen,intervallstellen);
      %      intervallwerte=intervallwerte(p,:);
      %  end
        
        A=zeros(m+ordnung(p),ordnung(p)+m);
        for j=1:m+ordnung(p)
            for q=1:ordnung(p)
                A(j,q)=(intervallstellen(j)-x1((i)+1))^(q-1)/faktorielle(q);
            end
            
            
            for sp=1:m
                
                if ordp>0
                    %A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*polyval(psi(sp,:,ordnung(p)),(intervallstellen(j)-x1(i+1))/h(i+1));
                    A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*psival(ordnung(p),sp,j);
                else                    
                    A(j,sp)=polyval(psi0(sp,:),(intervallstellen(j)-x1((i)+1))/h((i)+1));
                end
            end
        end
        help=A\intervallwerte';
        intervallstart((i)+1,1:length(help),p)=help;
    end
end

j=0;
start=zeros(N*max(ordnung)*n+N*m*n+parameter,1);
for i=0:N-1
    for k=1:max(ordnung)
        for p=1:n
            if k<=ordnung(p)
                j=j+1;start(j,1)=intervallstart((i)+1,k,p);
            end
        end
    end
    for k=1:m        
        for p=1:n
           j=j+1;start(j,1)=intervallstart((i)+1,k+ordnung(p),p);
        end
    end
end
for pii=1:parameter
    j=j+1;start(j,1)=par(pii);
end
start=start(1:j,1);
