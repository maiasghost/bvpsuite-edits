function output = equations(was,a,bvpfile,x1,wert,praediktor,psival,psi,linear)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equations prepares the equations and differential matrix for the Newton solver;
% % % FUNCTION CALL:
%output=equations(what,a,bvpfile,x1,value,praediktor,psival,psi,linear)
% INPUTS:  what          ... Determines what you would like to do with this
%                            routine
%                            'F' ... internal function call by run to
%                            determine the equations that will be solved by
%                            the Newton Solver
%                            'DF' ... the matrix that determines the
%                            direction for the Newton Solver (Differential
%                            Matrix of F)
%                            'residual' ... is called by the mesh
%                            adaptation and returns the residual
%                            'polynom' ... output of polynomial
%                            coefficients in monomial basis
%                            'valx1' ... calculation of the solution values on the
%                            vector x1 (mesh points)
%                            'valx1tau' calculation of the solution values
%                            on the vector x1 (mesh points) and tau
%                            (collocation points)
%                            'x1tau' ... output of mesh points and collocation
%                            points
%                            'basispolynome0' ... output of Runge Kutta basis
%                            coefficients for order 0
%                            'basispolynome' ... output of Runge Kutta
%                            basis coefficitens for other orders than 0
%                            'rho' ... output of the position of the
%                            collocation points transformed to the interval
%                            [0 1]
%                            'wert' ... returns special values of the
%                            solution on a defined vector
%                            'wertabl' ... the same for the first
%                            derivative
%                            'tau' ... returns the values of the collocation points 
%          a            ...  vector on which all the operations are done
%                           (when 'F' or 'DF' are chosen it is the vector on which the
%                            actual Newton step will be calculated)
%          wert        ...  the special vector for the 'what' options value
%                            and diffvalue
%          x1           ...  mesh of initial profile
%          praediktor   ...  for pathfollowing this variable gets the start vector in the first row
%                            the tangential vectori n the second row. If
%                            there is an additional entry in 
%                            praediktor(3,1) then the
%                            sum(ordnung)+1  equation is replaced by p1=praediktor(3,1)
%          psival,psi   ...  both stay the same for all the Newton
%                            iterations; they can be left out, but would
%                            then be calculated for each step which would
%                            be much slower
%          linear       ...  boolean; determines if problem is solved by
%                            linear solver or not
% % OUTPUTS: output       ...  totally depending on 'what' you want, normally
%                            a vector
% AUTHOR: Georg Kitzhofer
% ADAPTED BY: Christa Simon
% DATE: ??
% COMMENT: see the bvpsuite manual for the mathematical background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<9) linear=0; if (nargin<8) psi=[]; if (nargin<7) psival=[]; if (nargin<6) praediktor=[];  end;end;end;end;


%bf is a short form for bvpfile
bf=bvpfile;
%----------------------------------------------
%Call data from bvpfile

N=length(x1)-1;
b=x1((N)+1);
a1=x1((0)+1);
ordnung=feval(bf,'ordnung');
n=length(ordnung);%number of equations
rho=feval(bf,'rho');%position of the collocation points transformed to the interval [0 1]
standard=feval(bf,'standard');%are standard collocation points chosen? (Gauss, Lobatto, equidistant)
parameter=feval(bf,'parameter');%additional parameters
c=feval(bf,'c');%boundary values inside the interval
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];
%Determines the intervals where the points for the Additional Conditions
%are located
for i=1:length(c)
   j=0;
   if (c(i)<a1) || (c(i)>b)
       err('equations_err1');
   end
   while c(i)>=x1((j)+1) && j<N
      j=j+1;
   end
   cint(i)=j-1;
end
           if (standard)
   if (rho(1)==1)
       [rho]=gauss(rho(2));
   end
   if (rho(1)==2)
       [rho]=lobatto(rho(2));
   end
end

m=length(rho);%number of collocation points;

%----------------------------------------------
% Definition of h and the collocation points tau

h(1:N)=x1(2:N+1)-x1(1:N);
%in the rows of x1(1:N).'*ones(1,m) is xi
% in the columns are xi, ...xN
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
%rows of tau reflect the collocation points in i-th interval
%-----------------------------------------------
%Calculation of Runge-Kutta Basis
%
%
%psi(:,:,i) for n-th order psi(:,:,n-i) equals i-th derivative
%for example for 3-rd order:
%psi(:,:,1) ... second derivative
%psi(:,:,2) ... first derivative
%psi(:,:,3) ... no derivative

% psi contains the coefficients of the basisvectors
% polyval evaluates psi at the positions rho,1 and 0
% m = length(rho) ... number of collocation points

if length(psi)==0% is psi defined by user?
   for ord=1:max(ordnung)
       for i=1:m
           psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=Psi(i,rho,ord);
       end
   end
end
if length(psival)==0%is psi defined by user?
   for ord=1:max(ordnung)
       for i=1:m
           %evaluation of psi
           psival(ord,i,1:m)=polyval(psi(i,:,ord),rho(1:m));
           psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
           psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
       end
   end
end
   
%For 0-th order special definitions for psi are necessary

if min(ordnung)==0 || (length(was)==7 && length(strfind(was,'wertabl'))>0) || (length(was)==8 && length(strfind(was,'residual'))>0)
   for i=1:m
       psi0(i,:,1)=Psi(i,rho,0);
   end
end



%-----------------------------------------------
%entries of "a" will be transformed to the names of the documentation
%y,z,p see manual!
if ~linear
   j=0;y=0;z=0;p=0;
   for i=0:N-1
       for q=1:max(ordnung)
           for ni=1:n
               if q<=ordnung(ni)
                   j=j+1;y(ni,q,(i)+1)=a(j);
               end
           end
       end
       z(1:n,1:m,(i)+1)=reshape(a(j+1:j+n*m),n,m);
       j=j+n*m;
   end
   p(1:parameter)=a(j+1:j+parameter).';
else
   y=zeros(n,max(ordnung),N);
   z=zeros(n,m,N);
   p=zeros(parameter);
end


%-----------------------------------------------
% Calculation of boundary values
% if c is given the additional conditions are used
% all names of the variables are conform to the documentation
% (x1,h,y,z,psi,P)
% Poptimiert ... a faster evalutaion of P (for P itsself see documentation)

       switch was
  case 'F'

if ~linear
   output=zeros(length(a),1);
   j=0;
   if isempty(c)
       fza=0;fzb=0;
       for oi=1:max(ordnung)
           for ni=1:n
               if oi<=ordnung(ni)
                   fza(ni,oi)=y(ni,oi,(0)+1);
                   fzb(ni,oi)=Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
               end
           end
       end
       Rx=feval(bf,'R',[],[],[],[],fza,fzb,[],[],p);
       if sum(ordnung)+parameter>0
           output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
           j=j+sum(ordnung)+parameter;
       end
   else
       fzc=zeros(n,max(ordnung),length(c));
       for oi=1:max(ordnung)
           for ni=1:n
               if oi<=ordnung(ni)
                   for ci=1:length(c)
                       fzc(ni,oi,ci)=P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                   end
               end
           end
       end
       Rx=feval(bf,'R',[],[],[],[],[],[],[],[],p,[],[],[],[],fzc);

       if sum(ordnung)+parameter>0
           output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
           j=j+sum(ordnung)+parameter;
       end
   end


   if length(praediktor)>0
       if length(praediktor(:,1))==2
           output(sum(ordnung)+1,1)=sum((a.'-praediktor(1,:)).*praediktor(2,:));
       end
       if length(praediktor(:,1))==3
           output(sum(ordnung)+1,1)=p(1)-praediktor(3,1);
       end
   end

    %-----------------------------------------------
   %Evaluation of the equations for the solver
   fz=zeros(n,max(ordnung)+1);
   for i=0:N-1
       for k=1:m
           for oi=0:max(ordnung)-1
               for ni=1:n
                   if oi<ordnung(ni)
                       fz(ni,oi+1)=Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k); % psi wird immer von der Ordnung der entsprechenden Gleichung gewählt
                   end
               end
           end
           for ni=1:n
               fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
           end
           if n>0
               g_=feval(bf,'g',[],[],tau((i)+1,k),fz,[],[],[],[],p);
               output(j+1:j+n,1)=g_(1:n);
               j=j+n;
           end
       end
       if (i<N-1)
           for oi=0:max(ordnung)-1
               for ni=1:n
                   if oi<=ordnung(ni)-1
                       j=j+1;u=Poptimiert(oi,i,ni,x1((i+1)+1),m,x1,h,y,z,psival,ordnung(ni),m+1)-y(ni,oi+1,(i+1)+1);
                       output(j,1)=u;
                   end
               end
           end
       end
   end
else
   %Linear----------------------------------------------------------------
   output=zeros(length(a),1);
   j=0;
   if isempty(c)
       fza=zeros(n,max(ordnung));fzb=zeros(n,max(ordnung));
       Rx=feval(bf,'R',[],[],[],[],fza,fzb,[],[],p);
       if sum(ordnung)+parameter>0
           output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
           j=j+sum(ordnung)+parameter;
       end
   else
       fzc=zeros(n,max(ordnung),length(c));
       Rx=feval(bf,'R',[],[],[],[],[],[],[],[],p,[],[],[],[],fzc);
       if sum(ordnung)+parameter>0
           output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
           j=j+sum(ordnung)+parameter;
       end
   end
   if length(praediktor)>0
       if length(praediktor(:,1))==2
           output(sum(ordnung)+1,1)=sum((a.'-praediktor(1,:)).*praediktor(2,:));
       end
       if length(praediktor(:,1))==3
           output(sum(ordnung)+1,1)=p(1)-praediktor(3,1);
       end
   end

   %-----------------------------------------------
   %Evaluation of the equations for the solver
   fz=zeros(n,max(ordnung)+1);
   for i=0:N-1
       for k=1:m
           if n>0
               g_=feval(bf,'g',[],[],tau((i)+1,k),fz,[],[],[],[],p);
               output(j+1:j+n,1)=g_(1:n);
               j=j+n;
           end
       end
       if (i<N-1)
           for oi=0:max(ordnung)-1
               for ni=1:n
                   if oi<=ordnung(ni)-1
                       j=j+1;
                       output(j,1)=0;
                   end
               end
           end
       end
   end
   %End linear-------------------------------------------------------------
end
   %Calculate residual for defined vector
      case 'residual'
   respoints=wert;
   fz=zeros(n,max(ordnung)+1);
        for index=1:length(respoints)
            %Determine the interval where the residualpoint is located
            if (respoints(index)<x1(1)) || (respoints(index)>x1((N)+1))
                err('equations_err2');
            end
            for j=0:N
                if respoints(index)>=x1(j+1)
                    i=j;
                end
                if i==N
                    i=i-1;
                end
            end
           %Determine values of polynomials of the resiudalpoint in the
           %corresponding interval
           for oi=0:max(ordnung)
               for ni=1:n
                   if oi<ordnung(ni)+1
                       if oi<ordnung(ni)
                           fz(ni,oi+1)=P(oi,i,ni,respoints(index),m,x1,h,y,z,psi,ordnung(ni));
                       else
                           fz(ni,oi+1)=Pablspezial(oi,i,ni,respoints(index),m,x1,h,y,z,psi0,ordnung(ni));
                       end
                   end
               end
           end
           g_=feval(bf,'g',[],[],respoints(index),fz,[],[],[],[],p);
               output(1:n,index)=g_(1:n);
       end


%-----------------------------------------------


   case 'DF'

%JACOBIAN

if ~linear
%----Start Boundary conditions
%R_0
j=0;
% Dimensioning of non-zero entries
% Boundary conditions + Matrix
nonzeroent=(sum(ordnung)+parameter)*(max(ordnung)*n*max(2,length(c))+m*n+parameter)+...
          N*m*n*max(ordnung)*n+...
          N*m*n*m*n+...
          N*max(ordnung)*n*max(ordnung)*n+...
          N*max(ordnung)*n*m+...
          N*max(ordnung)*n*max(ordnung)*n;
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
OUTPUT=zeros(nonzeroent,1);
gen_c=0; %General count for sparse matrix efficiency
if isempty(c)
   fza=0;
   fzb=0;
   for oi=1:max(ordnung)
       for ni=1:n
           if oi<=ordnung(ni)
               fza(ni,oi)=y(ni,oi,(0)+1);
               fzb(ni,oi)=Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
           end
       end
   end
    % p1 index for boundary conditions of the equations, oi index for the
   % derivatives
   %p2 the "respect variables" of the differentiation have more than one
   %component
   DpR=feval(bf,'DpR',[],[],[],[],fza,fzb,[],[],p,[],[]);
   for p1=1:sum(ordnung)+parameter
       for k=1:max(ordnung)
           DRa_(k,:)=feval(bf,'DR',strcat('a',num2str(k)),p1,[],[],fza,fzb,[],[],p);
           DRb_(k,:)=feval(bf,'DR',strcat('b',num2str(k)),p1,[],[],fza,fzb,[],[],p);
       end
       if (N~=1)
           j=j+1;
           zaehler=0;
           zaehler2=0;
           %the derivatives with respect to y_0k follow
           for k=1:max(ordnung)
               help=feval(bf,'DR',strcat('a',num2str(k)),p1,[],[],fza,fzb,[],[],p);
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler=zaehler+1;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=zaehler;
                       OUTPUT(gen_c)=help(p2);
                   end
               end
           end
           %the derivatives with respect to y_(N-1)k follow
           for k=1:max(ordnung)
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler2=zaehler2+1;
                       help=0;
                       for abl=0:k-1
                           help=help+DRb_(abl+1,p2)*((b-x1((N-1)+1))^(k-abl-1)/faktorielle(k-abl));
                       end
                       if (k-abl-1)>=0
                           col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=col;
                           OUTPUT(gen_c)=help;
                       else
                           col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=col;
                           OUTPUT(gen_c)=0;
                       end
                   end
               end
           end
       else
           j=j+1;
           zaehler2=0;
           for k=1:max(ordnung)
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler2=zaehler2+1;
                       help=0;
                       for abl=0:k-1
                           help=help+DRb_(abl+1,p2)*((b-a1)^(k-abl-1)/faktorielle(k-abl));
                       end
                       if (k-abl-1)>=0
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=zaehler2;
                           OUTPUT(gen_c)=DRa_(k,p2)+help;
                       else
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=zaehler2;
                           OUTPUT(gen_c)=DRa_(k,p2);
                       end
                   end
               end
           end
       end
       %derivatives with respect to z_(N-1)q
       for q=1:m
           for p2=1:n
               zaehler2=zaehler2+1;
               help=0;
               for abl=0:ordnung(p2)-1
                   help=help+DRb_(abl+1,p2)*h((N-1)+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),1);
               end
               col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
               gen_c=gen_c+1;
               ROW(gen_c)=j;
               COL(gen_c)=col;
               OUTPUT(gen_c)=help;
           end
       end
       %derivatives with respect to the unknown parameters p_1 to
       %p_parameter
       if parameter>0
           col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
           gen_c_start=gen_c+1;
           gen_c=gen_c+length(col);
           ROW(gen_c_start:gen_c)=j;
           COL(gen_c_start:gen_c)=col;
           OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
       end
   end
else   % for additional conditions
   fzc=0;
   for oi=1:max(ordnung)
       for ni=1:n
           if oi<=ordnung(ni)
               for ci=1:length(c)
                   fzc(ni,oi,ci)=P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
               end
           end
       end
   end
   % p1 index for boundary conditions of the equations, oi index for the
   % derivatives
   %p2 the "respect variables" of the differentiation have more than one
   %component
   DpR=feval(bf,'DpR',[],[],[],[],[],[],[],[],p,[],[],[],[],fzc);
   for p1=1:sum(ordnung)+parameter
       for k=1:max(ordnung)
           for ci=1:length(c)
               DR_(k,:,ci)=feval(bf,'DR',strcat('c',num2str(ci),'_',num2str(k)),p1,[],[],[],[],[],[],p,[],[],[],[],fzc);
           end
       end
         % Derivation with respect to y_c1_k, z_c1_k, y_c2_k, z_c2_k, ...
          % -----------------------------------------------------------
       j=j+1;
       for ci=1:length(c)
           zaehler2=0;
           for k=1:max(ordnung)
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler2=zaehler2+1;
                       help2=0;
                        % Every polynomial containing the variable y_ci has
                       % to be differentiated and added
                       for cii=1:length(c)
                           if cint(cii)==cint(ci)
                               help=0;
                               for abl=0:k-1
                                   help=help+DR_(abl+1,p2,cii)*((c(cii)-x1((  cint(cii)  )+1))^(k-abl-1)/faktorielle(k-abl));
                               end
                               help2=help2+help;
                           end
                       end
                       col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=col;
                       OUTPUT(gen_c)=help2;
                   end
               end
           end
           %Derivation with respect to z_(N-1)q
           for q=1:m
               for p2=1:n
                   zaehler2=zaehler2+1;
                   help2=0;
                   for cii=1:length(c)
                       if cint(cii)==cint(ci)
                           help=0;
                           for abl=0:ordnung(p2)-1
                               help=help+DR_(abl+1,p2,cii)*h((cint(cii))+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),(c(cii)-x1((  cint(cii)  )+1))/h( cint(cii) +1));
                           end
                           help2=help2+help;
                       end
                   end
                   col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                   gen_c=gen_c+1;
                   ROW(gen_c)=j;
                   COL(gen_c)=col;
                   OUTPUT(gen_c)=help2;
               end
           end
       end

%-----------------------------------------------------------                     %Derivatives with respect to unknown parameters p_1 bis p_parameter
       if parameter>0
           col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
           gen_c_start=gen_c+1;
           gen_c=gen_c+length(col);
           ROW(gen_c_start:gen_c)=j;
           COL(gen_c_start:gen_c)=col;
           OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
       end
   end
end

%----------End Boundary Conditions

%----------Start J_i with C_i (notation according to manual)
fz=zeros(n,max(ordnung)+1);
for i=0:N-1
   for k=1:m
       for oi=0:max(ordnung)-1
           for ni=1:n
               if oi<ordnung(ni)
                   fz(ni,oi+1)=Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k);
               end
           end
       end
       for ni=1:n
           fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
       end
       %p1 indicates, which component of g_ik will be differentiated with respect
       %to y und z
       for oi=1:max(ordnung)+1
           Dg_(:,:,oi)=feval(bf,'Dg',oi,[],tau((i)+1,k),fz,[],[],[],[],p);
       end
       Dpg=feval(bf,'Dpg',[],[],tau((i)+1,k),fz,[],[],[],[],p);
       for p1=1:n
          j=j+1;
          %Differentiation with respect to y_i, sequence: y_i1 to
          %y_i(max(ordnung))
          zaehler=0;
          for oi2=1:max(ordnung)
             %p2 indicates, which component of y_i(oi2) will be
             %differentiated
             for p2=1:n
                 %Differentiation with respect to y_i(oi2) only if
                 %maximum order of the equation equals oi2
                 if oi2<=ordnung(p2)
                     zaehler=zaehler+1;
                     help=0;
                     for r1=1:oi2
                         help=help+Dg_(p1,p2,r1)*(((tau((i)+1,k)-x1((i)+1))^(oi2-r1))/(faktorielle(oi2-r1+1)));
                     end
                     col=i*m*n+sum(ordnung)*i+zaehler;
                     gen_c=gen_c+1;
                     ROW(gen_c)=j;
                     COL(gen_c)=col;
                     OUTPUT(gen_c)=help;
                 end
             end
          end
          zaehler=zaehler+1;
          for p2=1:n
              help(p2,1:m)=0;
              for r1=1:ordnung(p2)
                  help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,r1)*h((i)+1)^(ordnung(p2)-r1+1).*psival(ordnung(p2)-r1+1,1:m,k);
              end
              help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,ordnung(p2)+1).*(k==(1:m));
          end
          help2=help(:).';
          col=i*m*n+sum(ordnung)*i+zaehler:i*m*n+sum(ordnung)*i+zaehler+n*m-1;
          gen_c_start=gen_c+1;
          gen_c=gen_c+length(col);
          ROW(gen_c_start:gen_c)=j;
          COL(gen_c_start:gen_c)=col;
          OUTPUT(gen_c_start:gen_c)=help2;
          zaehler=zaehler+n*m;
          if parameter>0
              col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
              gen_c_start=gen_c+1;
              gen_c=gen_c+length(col);
              ROW(gen_c_start:gen_c)=j;
              COL(gen_c_start:gen_c)=col;
              OUTPUT(gen_c_start:gen_c)=Dpg(p1,1:parameter);
          end
       end
   end
      if i<N-1
       for abl=0:max(ordnung)-1
           for p1=1:n
               if abl<=ordnung(p1)-1
                   j=j+1;
                   %Differentiation with respect to y_i(oi2)
                   zaehler=0;
                   for oi2=1:max(ordnung)
                       for p2=1:n
                           if oi2<=ordnung(p2)
                               zaehler=zaehler+1;
                               if p1==p2
                                   if oi2-abl-1>=0
                                       col=i*(m*n+sum(ordnung))+zaehler;
                                       gen_c=gen_c+1;
                                       ROW(gen_c)=j;
                                       COL(gen_c)=col;
                                       OUTPUT(gen_c)=((x1((i+1)+1)-x1((i)+1))^(oi2-abl-1))/(faktorielle(oi2-abl));
                                   else
                                       col=i*(m*n+sum(ordnung))+zaehler;
                                       gen_c=gen_c+1;
                                       ROW(gen_c)=j;
                                       COL(gen_c)=col;
                                       OUTPUT(gen_c)=0;                                                                           end
                               end
                           end
                       end
                   end
                   %D with respect to z_iq
                   zaehler=0;
                   for q=1:m
                       zaehler=zaehler+1;
                       col=i*m*n+(i+1)*sum(ordnung)+(zaehler-1)*n+p1;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=col;
                       OUTPUT(gen_c)=h((i)+1)^(ordnung(p1)-abl)*psival(ordnung(p1)-(abl+1)+1,q,m+1);
                   end
                   zaehler=0;
                   for oi2=1:max(ordnung)
                       for p2=1:n
                           if oi2<=ordnung(p2)
                               zaehler=zaehler+1;
                               if p1==p2
                                   col=(i+1)*m*n+(i+1)*sum(ordnung)+zaehler;
                                   gen_c=gen_c+1;
                                   ROW(gen_c)=j;
                                   COL(gen_c)=col;
                                   OUTPUT(gen_c)=-(oi2==abl+1);
                               end
                           end
                       end
                   end
               end
           end
       end
   end
end

else
%Start linear-----------------
%----Start Boundary Conditions
%R_0
j=0;
%Dimensioning of non-zero-entries
%Boundary Conditions + Matrix
nonzeroent=(sum(ordnung)+parameter)*(max(ordnung)*n*max(2,length(c))+m*n+parameter)+...
   N*m*n*max(ordnung)*n+...
   N*m*n*m*n+...
   N*max(ordnung)*n*max(ordnung)*n+...
   N*max(ordnung)*n*m+...
   N*max(ordnung)*n*max(ordnung)*n;
ROW=zeros(nonzeroent,1);
COL=zeros(nonzeroent,1);
OUTPUT=zeros(nonzeroent,1);
gen_c=0; %General count for sparse matrix efficiency
if isempty(c)
   fza=zeros(n,max(ordnung));
   fzb=zeros(n,max(ordnung));
   %Used Notations: p1: indicates the boundary conditions of the equations, oi: indicates the
   %orders of the derivatives
   %p2: indicates the components of the "respect-varibles"
   DpR=feval(bf,'DpR',[],[],[],[],fza,fzb,[],[],p,[],[]);
   for p1=1:sum(ordnung)+parameter
       for k=1:max(ordnung)
           DRa_(k,:)=feval(bf,'DR',strcat('a',num2str(k)),p1,[],[],fza,fzb,[],[],p);
           DRb_(k,:)=feval(bf,'DR',strcat('b',num2str(k)),p1,[],[],fza,fzb,[],[],p);
       end
       if (N~=1)
           j=j+1;
           zaehler=0;
           zaehler2=0;
           %derivatives with respect to y_0k
           for k=1:max(ordnung)
               help=feval(bf,'DR',strcat('a',num2str(k)),p1,[],[],fza,fzb,[],[],p);
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler=zaehler+1;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=zaehler;
                       OUTPUT(gen_c)=help(p2);
                   end
               end
           end
          %derivatives with respect to y_(N-1)k
           for k=1:max(ordnung)
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler2=zaehler2+1;
                       help=0;
                       for abl=0:k-1
                           help=help+DRb_(abl+1,p2)*((b-x1((N-1)+1))^(k-abl-1)/faktorielle(k-abl));
                       end
                       if (k-abl-1)>=0
                           col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=col;
                           OUTPUT(gen_c)=help;
                       else
                           col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=col;
                           OUTPUT(gen_c)=0;
                       end
                   end
               end
           end
       else
           j=j+1;
           zaehler2=0;
           for k=1:max(ordnung)
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler2=zaehler2+1;
                       help=0;
                       for abl=0:k-1
                           help=help+DRb_(abl+1,p2)*((b-a1)^(k-abl-1)/faktorielle(k-abl));
                       end
                       if (k-abl-1)>=0
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=zaehler2;
                           OUTPUT(gen_c)=DRa_(k,p2)+help;
                       else
                           gen_c=gen_c+1;
                           ROW(gen_c)=j;
                           COL(gen_c)=zaehler2;
                           OUTPUT(gen_c)=DRa_(k,p2);
                       end
                   end
               end
           end
       end
       %derivatives with respect to z_(N-1)q
       for q=1:m
           for p2=1:n
               zaehler2=zaehler2+1;
               help=0;
               for abl=0:ordnung(p2)-1
                   help=help+DRb_(abl+1,p2)*h((N-1)+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),1);
               end
               col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
               gen_c=gen_c+1;
               ROW(gen_c)=j;
               COL(gen_c)=col;
               OUTPUT(gen_c)=help;
           end
       end
       %derivatives with respect to the parameteres p_1 to p_parameter
       if parameter>0
           col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
           gen_c_start=gen_c+1;
           gen_c=gen_c+length(col);
           ROW(gen_c_start:gen_c)=j;
           COL(gen_c_start:gen_c)=col;
           OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
       end
   end
else   % The same for additional conditions
   fzc=zeros(n,max(ordnung),length(c));
   %Used Notations: p1: indicates the boundary conditions of the equations, oi: indicates the
   %orders of the derivatives
   %p2: indicates the components of the "respect-varibles"
   DpR=feval(bf,'DpR',[],[],[],[],[],[],[],[],p,[],[],[],[],fzc);
   for p1=1:sum(ordnung)+parameter
       for k=1:max(ordnung)
           for ci=1:length(c)
               DR_(k,:,ci)=feval(bf,'DR',strcat('c',num2str(ci),'_',num2str(k)),p1,[],[],[],[],[],[],p,[],[],[],[],fzc);
           end
       end
       % derivations with respect to y_c1_k, z_c1_k, y_c2_k, z_c2_k
       % -----------------------------------------------------------
       j=j+1;
       for ci=1:length(c)
           zaehler2=0;
           for k=1:max(ordnung)
               for p2=1:n
                   if k<=ordnung(p2)
                       zaehler2=zaehler2+1;
                       help2=0;
                       %every polynomial containing y_ci has to be
                       %differentiated
                       for cii=1:length(c)
                           if cint(cii)==cint(ci)
                               help=0;
                               for abl=0:k-1
                                   help=help+DR_(abl+1,p2,cii)*((c(cii)-x1((  cint(cii)  )+1))^(k-abl-1)/faktorielle(k-abl));
                               end
                               help2=help2+help;
                           end
                       end
                       col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=col;
                       OUTPUT(gen_c)=help2;
                   end
               end
           end
           %derivatives with respect to z_(N-1)q
           for q=1:m
               for p2=1:n
                   zaehler2=zaehler2+1;
                   help2=0;
                   for cii=1:length(c)
                       if cint(cii)==cint(ci)
                           help=0;
                           for abl=0:ordnung(p2)-1
                               help=help+DR_(abl+1,p2,cii)*h((cint(cii))+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),(c(cii)-x1((  cint(cii)  )+1))/h( cint(cii) +1));
                           end
                           help2=help2+help;
                       end
                   end
                   col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                   gen_c=gen_c+1;
                   ROW(gen_c)=j;
                   COL(gen_c)=col;
                   OUTPUT(gen_c)=help2;
               end
           end
       end

       %-----------------------------------------------------------

       %derivatives with respect to the parameters p_1 to p_parameter
       if parameter>0
           col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
           gen_c_start=gen_c+1;
           gen_c=gen_c+length(col);
           ROW(gen_c_start:gen_c)=j;
           COL(gen_c_start:gen_c)=col;
           OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
       end
   end
end

%----------End Boundary Conditions

%----------Start J_i with C_i (notation equals manual)
fz=zeros(n,max(ordnung)+1);
for i=0:N-1
   for k=1:m
       %p1 indicates, which component of g_ik will be differentiated with respect
       %to y und z
       for oi=1:max(ordnung)+1
           Dg_(:,:,oi)=feval(bf,'Dg',oi,[],tau((i)+1,k),fz,[],[],[],[],p);
       end
       Dpg=feval(bf,'Dpg',[],[],tau((i)+1,k),fz,[],[],[],[],p);
       for p1=1:n
           j=j+1;
           %Differentiation with respect to y_i, sequence: y_i1 to
           %y_i(max(ordnung))
           zaehler=0;
           for oi2=1:max(ordnung)
               %p2 indicates, which component of y_i(oi2) will be
               %differentiated
               for p2=1:n
                   %Differentiation with respect to y_i(oi2) only if
                   %maximum order of the equation equals oi2
                   if oi2<=ordnung(p2)
                       zaehler=zaehler+1;
                       help=0;
                       for r1=1:oi2
                           help=help+Dg_(p1,p2,r1)*(((tau((i)+1,k)-x1((i)+1))^(oi2-r1))/(faktorielle(oi2-r1+1)));
                       end
                       col=i*m*n+sum(ordnung)*i+zaehler;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=col;
                       OUTPUT(gen_c)=help;
                   end
               end
           end
           zaehler=zaehler+1;
           for p2=1:n
               help(p2,1:m)=0;
               for r1=1:ordnung(p2)
                   help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,r1)*h((i)+1)^(ordnung(p2)-r1+1).*psival(ordnung(p2)-r1+1,1:m,k);
               end
               help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,ordnung(p2)+1).*(k==(1:m));
           end
           help2=help(:).';
           col=i*m*n+sum(ordnung)*i+zaehler:i*m*n+sum(ordnung)*i+zaehler+n*m-1;
           gen_c_start=gen_c+1;
           gen_c=gen_c+length(col);
           ROW(gen_c_start:gen_c)=j;
           COL(gen_c_start:gen_c)=col;
           OUTPUT(gen_c_start:gen_c)=help2;
           zaehler=zaehler+n*m;
           if parameter>0
               col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
               gen_c_start=gen_c+1;
               gen_c=gen_c+length(col);
               ROW(gen_c_start:gen_c)=j;
               COL(gen_c_start:gen_c)=col;
               OUTPUT(gen_c_start:gen_c)=Dpg(p1,1:parameter);
           end
       end
   end

   if i<N-1
       for abl=0:max(ordnung)-1
           for p1=1:n
               if abl<=ordnung(p1)-1
                   j=j+1;
                   %derviatives with respect to y_i(oi2)
                   zaehler=0;
                   for oi2=1:max(ordnung)
                       for p2=1:n
                           if oi2<=ordnung(p2)
                               zaehler=zaehler+1;
                               if p1==p2
                                   if oi2-abl-1>=0
                                       col=i*(m*n+sum(ordnung))+zaehler;
                                       gen_c=gen_c+1;
                                       ROW(gen_c)=j;
                                       COL(gen_c)=col;
                                       OUTPUT(gen_c)=((x1((i+1)+1)-x1((i)+1))^(oi2-abl-1))/(faktorielle(oi2-abl));
                                   else
                                       col=i*(m*n+sum(ordnung))+zaehler;
                                       gen_c=gen_c+1;
                                       ROW(gen_c)=j;
                                       COL(gen_c)=col;
                                       OUTPUT(gen_c)=0;
                                   end
                               end
                           end
                       end
                   end
                   %D with respect to z_iq
                   zaehler=0;
                   for q=1:m
                       zaehler=zaehler+1;
                       col=i*m*n+(i+1)*sum(ordnung)+(zaehler-1)*n+p1;
                       gen_c=gen_c+1;
                       ROW(gen_c)=j;
                       COL(gen_c)=col;
                       OUTPUT(gen_c)=h((i)+1)^(ordnung(p1)-abl)*psival(ordnung(p1)-(abl+1)+1,q,m+1);
                   end
                   zaehler=0;
                   for oi2=1:max(ordnung)
                       for p2=1:n
                           if oi2<=ordnung(p2)
                               zaehler=zaehler+1;
                               if p1==p2
                                   col=(i+1)*m*n+(i+1)*sum(ordnung)+zaehler;
                                   gen_c=gen_c+1;
                                   ROW(gen_c)=j;
                                   COL(gen_c)=col;
                                   OUTPUT(gen_c)=-(oi2==abl+1);
                               end
                           end
                       end
                   end
               end
           end
       end
   end
end
%End linear-------------------
end











dimsparse=N*(n*m+sum(ordnung))+parameter;
ROW=ROW(1:gen_c);
COL=COL(1:gen_c);
OUTPUT=OUTPUT(1:gen_c);
output=sparse(ROW,COL,OUTPUT,dimsparse,dimsparse);

%--------End J_i with C_i

%Change of the derivation of the sum(ordnung)+1 th row when using
%pathfollowing

if length(praediktor)>0
   if length(praediktor(:,1))==2
       output(sum(ordnung)+1,:)=praediktor(2,:);
   end
   if length(praediktor(:,1))==3
       output(sum(ordnung)+1,N*(n*m+sum(ordnung))+1)=1;
   end
end



%--------Ende Jacobian
%Jacobi=output

  case 'polynom'

zeile=0;       for i=0:N-1
   for komp=1:n
       zeile=zeile+1;
       if ordnung(komp)~=0
           help=evalP0(i,komp,m,x1,h,y,z,psi,ordnung(komp));
       else
           help=evalP0(i,komp,m,x1,h,[],z,psi0,0);
       end
       output(zeile,m+max(ordnung)+1-length(help):m+max(ordnung))=help;
   end
end

%C: Values of the polynomials
  case 'valx1'
for i=1:n
   for j=0:N-1
       if ordnung(i)~=0
           help=P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));
       else
           help=P(0,j,i,x1((j)+1),m,x1,h,[],z,psi0,0);
       end
       output(i,j+1)=help;
   end
   if ordnung(i)~=0
       help=P(0,N-1,i,x1((N)+1),m,x1,h,y,z,psi,ordnung(i));
   else
       help=P(0,N-1,i,x1((N)+1),m,x1,h,[],z,psi0,0);
   end
   output(i,N+1)=help;
end

 case 'valx1tau'
output=zeros(n,N+N*m+1);
for i=1:n
   stelle=0;
   for j=0:N-1
       if ordnung(i)~=0
           if j~=0
               help=Poptimiert(0,j,i,x1((j)+1),m,x1,h,y,z,psival,ordnung(i),m+2);
           else
               help=P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));
           end
       else
           help=P(0,j,i,x1((j)+1),m,x1,h,[],z,psi0,0);
       end
       stelle=stelle+1;output(i,stelle)=help;
       for k=1:m
           if ordnung(i)~=0
               help=Poptimiert(0,j,i,tau((j)+1,k),m,x1,h,y,z,psival,ordnung(i),k);
           else
               help=P(0,j,i,tau((j)+1,k),m,x1,h,[],z,psi0,0);
           end
           stelle=stelle+1;output(i,stelle)=help;
       end
   end
   if ordnung(i)~=0
       help=P(0,N-1,i,x1((N)+1),m,x1,h,y,z,psi,ordnung(i));
   else
       help=P(0,N-1,i,x1((N)+1),m,x1,h,[],z,psi0,0);
   end
   stelle=stelle+1;
   output(i,stelle)=help;
end

 case 'x1tau'
stelle=0;
for j=0:N-1
     stelle=stelle+1;output(stelle)=x1((j)+1);
     for k=1:m
         stelle=stelle+1;output(stelle)=tau((j)+1,k);
     end
end
stelle=stelle+1;output(stelle)=x1((N)+1);

 case 'basispolynome0'
     output=psi0;      case 'basispolynome'
     output=psi;
      case 'rho'
     output=rho;
      case 'wert'
    for k=1:length(wert)         
        if (wert(k)<x1(1)) || (wert(k)>x1((N)+1))
            disp(wert);
            warning('Value outside range');
        end
        for i=0:N
            if wert(k)>=x1((i)+1)
                j=i;
            end
            if j==N
                j=j-1;
            end
        end
        for ni=1:n
            if ordnung(ni)~=0
                help=P(0,j,ni,wert(k),m,x1,h,y,z,psi,ordnung(ni));
            else
                help=P(0,j,ni,wert(k),m,x1,h,[],z,psi0,0);
            end
            output(ni,k)=help;
        end
    end

  case 'wertabl'
    for k=1:length(wert)
        if (wert(k)<x1(1)) || (wert(k)>x1((N)+1))
            err('equations_err3');
        end
        for i=0:N
            if wert(k)>=x1((i)+1)
                j=i;
            end
            if j==N
                j=j-1;
            end
        end
        for ni=1:n
            if ordnung(ni)~=1
                help=P(1,j,ni,wert(k),m,x1,h,y,z,psi,ordnung(ni));
            else
                help=Pablspezial(1,j,ni,wert(k),m,x1,h,[],z,psi0,1);
            end
            output(ni,k)=help;
        end
    end
       case 'tau'
    output=tau;

end %end of Switch

%--------Local functions

%C: help=P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));


function ret=P(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynomial of derivation abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
   if j-abl-1>=0
       sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
   end
end
if ord==0
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
   end
else
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
   end
end
ret=sum;

function ret=Pablspezial(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynom der Ableitung abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
   if j-abl-1>=0
       sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
   end
end
if ord-abl==0
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
   end
else
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
   end
end
ret=sum;

function ret=Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k)
%Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^([abl+1:ord]-abl-1))./(faktorielle([abl+1:ord]-abl)).*y(komp,abl+1:ord,(i)+1))+...
   sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));


function ret=evalP0(i,komp,m,x1,h,y,z,psi,ord)
if (nargin<9) ord=2; end;

     sum=0;
     for j=1:ord

        poly=[1 -x1((i)+1)];
        for k=2:(j-1)
                poly=conv(poly,[1 -x1((i)+1)]);
        end
        if (j-1)==0
            poly=[1];
        end
        zwisch(m+ord-j+1:m+ord)=poly/(factorial(j-1))*y(komp,j,(i)+1);
        sum=sum+zwisch;
    end
for j=1:m
     if ord==0
         help=polyvalpoly(psi(j,:,1),[1 -x1((i)+1)]/h((i)+1));
     else
         help=polyvalpoly(psi(j,:,ord-1+1),[1 -x1((i)+1)]/h((i)+1));
     end
     %Attetion: help has to be of length m+ord - correction for
     %mixed order
     help=help(length(help)-(m+ord-1):length(help));
     zwisch(1:m+ord)=h((i)+1)^ord*(z(komp,j,(i)+1)*help);
     sum=sum+zwisch;
end
ret=sum;

function out=polyvalpoly(psi,E)
%polynomial as value in polynomial receiving polynomial

PsiK(1,1)=1;t=1;
for k=2:length(psi)
   t=conv(t,E);
   for m=1:length(t)
     PsiK(k,m)=t(m);
   end
end
PsiE(((length(E)-1))*(length(psi)-1)+1)=PsiK(1,1)*psi(length(psi));
PsiEout=PsiE;
for k=2:length(psi)
   PsiE(((length(E)-1))*(length(psi)-1)+1-(length(E)-1)*(k-1):((length(E)-1))*(length(psi)-1)+1)=...
                      psi(length(psi)-k+1)*PsiK(k,1:(length(E)-1)*(k-1)+1);
   PsiEout=PsiEout+PsiE;
end

out=PsiEout;   function Psireturn=Psi(n,rho,nr)
 i=length(rho);
 prod=[1];
 for s=1:i
     if (s~=n)
         prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
     end
 end
 for s=1:(nr)
     prod=polyint(prod);
 end
 Psireturn=prod;
     %tabulars for collocation points

function [rho]=gauss(k)

%psi(:,:,i) for n-th order psi(:,:,n-i) equals i-th derivative
%for example for 3-rd order:
%psi(:,:,1) ... second derivative
%psi(:,:,2) ... first derivative
%psi(:,:,3) ... no derivative

switch k
   case 1
  rho=[0.5];
  psi(1,:,1)=[0 1 0];
  psi(1,:,2)=[0.5 0 0];
      case 2
       %C: rho: ci in Gaussian collocation, zeros of Legendre Polynomials shiftet to [0,1]
  %C: transformation from [-1,1] to [0,1]: y= 1/2x+1/2
  rho=[0.211324865405187117745425609748  0.788675134594812882254574390252];
  psi(1,:,1)=[0 -0.866025403784438646763723170755  1.36602540378443864676372317076 0];
  psi(2,:,1)=[0  0.866025403784438646763723170755   -0.366025403784438646763723170754 0];
  %0.th derivative
  %firs coll.point
  psi(1,:,2)=[-0.288675134594812882254574390252  0.683012701892219323381861585378 0 0];
  %second coll.point
  psi(2,:,2)=[0.288675134594812882254574390252  -0.183012701892219323381861585377 0 0];
     case 3
  rho=[1/2-1/10*15^(1/2) 1/2 1/2+1/10*15^(1/2)];
  psi(1,:,1)=[0 10/9 -5/3-1/6*15^(1/2) 5/6+1/6*15^(1/2) 0];
  psi(2,:,1)=[0 -20/9 10/3 -2/3 0];
  psi(3,:,1)=[0 10/9 -5/3+1/6*15^(1/2) 5/6-1/6*15^(1/2) 0];
  psi(1,:,2)=[5/18 -5/9-1/18*15^(1/2) 5/12+1/12*15^(1/2) 0 0];
  psi(2,:,2)=[-5/9 10/9 -1/3 0 0];
  psi(3,:,2)=[5/18 -5/9+1/18*15^(1/2) 5/12-1/12*15^(1/2) 0 0];
     case 4

  rho=[1/2-1/70*(525+70*30^(1/2))^(1/2) 1/2-1/70*(525-70*30^(1/2))^(1/2) 1/2+1/70*(525-70*30^(1/2))^(1/2) 1/2+1/70*(525+70*30^(1/2))^(1/2)];
  psi(1,:,1)=[0 -245/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/36*(105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             -7/24*(45+(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0];
  psi(2,:,1)=[0 245/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/36*(105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             -7/24*(-45+30^(1/2)-(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0];
  psi(3,:,1)=[0 -245/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/36*(-105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             7/24*(-45+30^(1/2)+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(-35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0];
  psi(4,:,1)=[0 245/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/36*(-105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             7/24*(45-(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(-35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0];
  psi(1,:,2)=[-49/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/144*(105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             -7/72*(45+(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/240*(35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0 0];
  psi(2,:,2)=[49/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/144*(105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             7/72*(45-30^(1/2)+(525-70*30^(1/2))^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) ...
             1/240*(35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0 0];
  psi(3,:,2)=[-49/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/144*(-105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             7/72*(-45+30^(1/2)+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             1/240*(-35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0 0];
  psi(4,:,2)=[49/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/144*(-105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             7/72*(45-(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/240*(-35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0 0];

 case 5

  rho=[1/2-1/42*(245+14*70^(1/2))^(1/2) 1/2-1/42*(245-14*70^(1/2))^(1/2) 1/2 1/2+1/42*(245-14*70^(1/2))^(1/2) 1/2+1/42*(245+14*70^(1/2))^(1/2)];
  psi(1,:,1)=[0 567/25/(35+2*70^(1/2))*70^(1/2) -27/40*(84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/20*(343+9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/280*(1911+77*(245+14*70^(1/2))^(1/2)+42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/280*(21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0];
  psi(2,:,1)=[0 567/25/(-35+2*70^(1/2))*70^(1/2) -27/40*(84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/20*(-343+2*70^(1/2)-9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             3/280*(-1911+42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/280*(21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0];
  psi(3,:,1)=[0 336/25 -168/5 1232/45 -112/15 8/15 0];
  psi(4,:,1)=[0 567/25/(-35+2*70^(1/2))*70^(1/2) 27/40*(-84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/20*(-343+2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/280*(1911-42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             3/280*(-21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0];
  psi(5,:,1)=[0 567/25/(35+2*70^(1/2))*70^(1/2) 27/40*(-84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/20*(-343+9*(245+14*70^(1/2))^(1/2)-2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/280*(-1911+77*(245+14*70^(1/2))^(1/2)-42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/280*(-21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0];
  psi(1,:,2)=[189/50/(35+2*70^(1/2))*70^(1/2) -27/200*(84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/80*(343+9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -1/280*(1911+77*(245+14*70^(1/2))^(1/2)+42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/560*(21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0 0];
  psi(2,:,2)=[189/50/(-35+2*70^(1/2))*70^(1/2) -27/200*(84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             3/80*(343-2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             1/280*(-1911+42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/560*(21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0 0];
  psi(3,:,2)=[56/25 -168/25 308/45 -112/45 4/15 0 0];
  psi(4,:,2)=[189/50/(-35+2*70^(1/2))*70^(1/2) 27/200*(-84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) -3/80*(-343+2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
            -1/280*(1911-42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
            3/560*(-21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0 0];
  psi(5,:,2)=[189/50/(35+2*70^(1/2))*70^(1/2) 27/200*(-84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/80*(343-9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))*70^(1/2)/(35+2*70^(1/2)) ...
             1/280*(-1911+77*(245+14*70^(1/2))^(1/2)-42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/560*(-21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0 0];
        case 6

  rho=[0.33765242898423986094e-1 0.16939530676686774317 0.38069040695840154568 0.61930959304159845432 0.83060469323313225683 0.96623475710157601391];
  psi(1,:,1)=[0 -8.1412617290086756177 28.978672208695686824 -40.408362140839295291 27.785390555068868256 -9.6944498478780709320 1.5656732001510719330 0];
  psi(2,:,1)=[0 24.533738740695531559 -83.334379226361945114 107.81105264377978693 -64.863346973441684340 16.973778445028729194 -.94046284317634892902 0];
  psi(3,:,1)=[0 -36.168340495472103373 113.68329746902200059 -130.85406519770904404 65.101388388983683720 -12.145253252968680079 .61693005543048870860 0];
  psi(4,:,1)=[0 36.168340495472103373 -103.32674550381061965 104.96268528468059170 -44.848707621074554040 7.6576120141334378913 -.37922770211461375460 0];
  psi(5,:,1)=[0 -24.533738740695531559 63.868053217811244240 -59.145237622403034740 23.711846151968643422 -3.9123422341959200109 .19180001403866795482 0];
  psi(6,:,1)=[0 8.1412617290086756177 -19.868898165356366883 17.633927032490995438 -6.8865705015049570236 1.1206548758805039360 -.54712724329265912892e-1 0];
  psi(1,:,2)=[-1.1630373898583822311 4.8297787014492811373 -8.0816724281678590582 6.9463476387672170640 -3.2314832826260236440 .78283660007553596650 0 0];
  psi(2,:,2)=[3.5048198200993616513 -13.889063204393657519 21.562210528755957386 -16.215836743360421085 5.6579261483429097313 -.47023142158817446451 0 0];
  psi(3,:,2)=[-5.1669057850674433390 18.947216244837000098 -26.170813039541808808 16.275347097245920930 -4.0484177509895600263 .30846502771524435430 0 0];
  psi(4,:,2)=[5.1669057850674433390 -17.221124250635103275 20.992537056936118340 -11.212176905268638510 2.5525373380444792971 -.18961385105730687730 0 0];
  psi(5,:,2)=[-3.5048198200993616513 10.644675536301874040 -11.829047524480606948 5.9279615379921608555 -1.3041140780653066703 .95900007019333977410e-1 0 0];
  psi(6,:,2)=[1.1630373898583822311 -3.3114830275593944805 3.5267854064981990876 -1.7216426253762392559 .37355162529350131200 -.27356362164632956446e-1 0 0];
  case 7
  rho=[.254460438286207377369051579761e-1, .129234407200302780068067613360, .297077424311301416546696793962, .500000000000000000000000000000, .702922575688698583453303206038, .870765592799697219931932386640, .974553956171379262263094842024];
 case 8
  rho=[.198550717512318841582195657153e-1, .101666761293186630204223031762, .237233795041835507091130475405, .408282678752175097530261928820, .591717321247824902469738071180, .762766204958164492908869524595, .898333238706813369795776968238, .980144928248768115841780434285];
 case 9
  rho=[.159198802461869550822118985482e-1, .819844463366821028502851059651e-1, .193314283649704801345648980329, .337873288298095535480730992678, .500000000000000000000000000000, .662126711701904464519269007322, .806685716350295198654351019671, .918015553663317897149714894035, .984080119753813044917788101452];
 case 10
  rho=[.130467357414141399610179939578e-1, .674683166555077446339516557883e-1, .160295215850487796882836317443, .283302302935376404600367028417, .425562830509184394557586999435, .574437169490815605442413000565, .716697697064623595399632971583, .839704784149512203117163682557, .932531683344492255366048344212, .986953264258585860038982006042];
 case 11
  rho=[.108856709269715035980309994386e-1, .564687001159523504624211153480e-1, .134923997212975337953291873984, .240451935396594092037137165271, .365228422023827513834234007300, .500000000000000000000000000000, .634771577976172486165765992700, .759548064603405907962862834729, .865076002787024662046708126016, .943531299884047649537578884652, .989114329073028496401969000561];
 case 12
  rho=[.921968287664037465472545492536e-2, .479413718147625716607670669405e-1, .115048662902847656481553083394, .206341022856691276351648790530, .316084250500909903123654231678, .437383295744265542263779315268, .562616704255734457736220684732, .683915749499090096876345768322, .793658977143308723648351209470, .884951337097152343518446916606, .952058628185237428339232933060, .990780317123359625345274545075];
 case 13
  rho=[.790847264070592526358527559645e-2, .412008003885110173967260817496e-1, .992109546333450436028967552086e-1, .178825330279829889678007696502, .275753624481776573561043573936, .384770842022432602967235939451, .500000000000000000000000000000, .615229157977567397032764060549, .724246375518223426438956426064, .821174669720170110321992303498, .900789045366654956397103244791, .958799199611488982603273918250, .992091527359294074736414724404];
 case 14
  rho=[.685809565159383057920136664797e-2, .357825581682132413318044303111e-1, .863993424651175034051026286748e-1, .156353547594157264925990098490, .242375681820922954017354640724, .340443815536055119782164087916, .445972525646328168966877674890, .554027474353671831033122325110, .659556184463944880217835912084, .757624318179077045982645359276, .843646452405842735074009901510, .913600657534882496594897371325, .964217441831786758668195569689, .993141904348406169420798633352];
 case 15
  rho=[.600374098975728575521714070669e-2, .313633037996470478461205261449e-1, .758967082947863918996758396129e-1, .137791134319914976291906972693, .214513913695730576231386631373, .302924326461218315051396314509, .399402953001282738849685848303, .500000000000000000000000000000, .600597046998717261150314151697, .697075673538781684948603685491, .785486086304269423768613368627, .862208865680085023708093027307, .924103291705213608100324160387, .968636696200352952153879473855, .993996259010242714244782859293];   end

function [rho]=lobatto(k)
if k<2 || k>15
   err('equations_err4');
end
switch k
   case 2
       rho=[0,1];
   case 3
       rho=[0., .500000000000000000000000000000, 1.];
   case 4
       rho=[0., .276393202250021030359082633127, .723606797749978969640917366873, 1.];
   case 5
       rho=[0., .172673164646011428100853771876, .500000000000000000000000000000, .827326835353988571899146228124, 1.];
   case 6
       rho=[0., .117472338035267653574498513019, .357384241759677451842924502980, .642615758240322548157075497020, .882527661964732346425501486981, 1.];
   case 7
       rho=[0., .84888051860716535063983893017e-1, .265575603264642893098114059046, .500000000000000000000000000000, .734424396735357106901885940954, .915111948139283464936016106983, 1.];
   case 8
       rho=[0., .641299257451966923312771193897e-1, .204149909283428848927744634301, .395350391048760565615671369827, .604649608951239434384328630173, .795850090716571151072255365699, .935870074254803307668722880610, 1.00000000000000000000000000000];
   case 9
       rho=[0., .501210022942699213438273777908e-1, .161406860244631123277057286454, .318441268086910920644623965646, .500000000000000000000000000000, .681558731913089079355376034354, .838593139755368876722942713546, .949878997705730078656172622209, 1.00000000000000000000000000000];
   case 10
       rho=[0., .402330459167705930855336695888e-1, .130613067447247462498446912570, .261037525094777752169412453634, .417360521166806487686890117021, .582639478833193512313109882979, .738962474905222247830587546366, .869386932552752537501553087430, .959766954083229406914466330411, 1.00000000000000000000000000000];
   case 11
       rho=[0., .329992847959704328338629319503e-1, .107758263168427790688791091946, .217382336501897496764518015261, .352120932206530304284044242220, .500000000000000000000000000000, .647879067793469695715955757780, .782617663498102503235481984739, .892241736831572209311208908054, .967000715204029567166137068050, 1.00000000000000000000000000000];
   case 12
       rho=[0., .275503638885588882962099308484e-1, .903603391779966608256792091415e-1, .183561923484069661168797572778, .300234529517325533867825104217, .431723533572536222567969072130, .568276466427463777432030927870, .699765470482674466132174895783, .816438076515930338831202427222, .909639660822003339174320790858, .972449636111441111703790069152, 1.00000000000000000000000000000];
   case 13
       rho=[0., .233450766789180440515472676223e-1, .768262176740638415670371964506e-1, .156905765459121286963620480217, .258545089454331899126531383182, .375356534946880003715663149813, .500000000000000000000000000000, .624643465053119996284336850187, .741454910545668100873468616818, .843094234540878713036379519783, .923173782325936158432962803549, .976654923321081955948452732378, 1.00000000000000000000000000000];
   case 14
       rho=[0., .200324773663695493224499189923e-1, .660994730848263744998898985459e-1, .135565700454336929707663799740, .224680298535676472341688647070, .328637993328643577478048298179, .441834065558148066170611645132, .558165934441851933829388354868, .671362006671356422521951701821, .775319701464323527658311352930, .864434299545663070292336200260, .933900526915173625500110101454, .979967522633630450677550081008, 1.00000000000000000000000000000];
   case 15
       rho=[0., .173770367480807136020743039652e-1, .574589778885118505872991842589e-1, .118240155024092399647940762012, .196873397265077144438235030682, .289680972643163759539051530631, .392323022318102880887160276864, .500000000000000000000000000000, .607676977681897119112839723136, .710319027356836240460948469369, .803126602734922855561764969318, .881759844975907600352059237988, .942541022111488149412700815741, .982622963251919286397925696035, 1.00000000000000000000000000000];
end