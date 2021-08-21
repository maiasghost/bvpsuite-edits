function [x1,start]=initialmesh(bvpfile,stellen,werte,x1,par,lambda,errest_call)
%INITIALMESH interpolates given values with cubic spline interpolation
%and calculates the corresponding coefficients in the Runge-Kutta basis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALL:
%[x1,start]=initialmesh(bvpfile,stellen,werte,x1,par)
%INPUT:
%bvpfile               stored bvpfile; contains the problem defining
%                      information
%stellen               mesh points on which the solution is known; entered
%                      by the user (=initial mesh)
%werte                 solution values corresponding to the given mesh
%                      points "stellen"
%x1                    new mesh; integer if semi-infinite interval
%par                   starting parameters
%lambda                starting eigenvalue (only relevant for Eigenvalue Problems)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT:that
%x1                    corresponding to the mesh x1
%                      Runge-Kutta coeff. are calculated
%start                 starting coefficients in the Runge-Kutta basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AUTHOR: gkitzhofer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADAPTED BY: csimon. 04/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DATE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMMENT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if (nargin<7)
koeff=[];
%
%       if (nargin<5)
%         par=[];
%     end;
% end;


if (nargin<7)
  errest_call=0;
  if (nargin<6)
    lambda =[];
    if (nargin<5)
      par=[];
    end
  end
end


N=length(x1)-1;

if errest_call
  
  n=feval(bvpfile,'n');
  
  ordnung=feval(bvpfile,'ordnung');
  
  
else
  
  if feval(bvpfile,'Infinite') &&  feval(bvpfile,'EVP')==0
    %if the problem is posed on a semi-finite interval but no eigenvalue
    %problem
    
    n=feval(bvpfile,'n'); %dimension of the system of ODEs
    k=n;
    
    if feval(bvpfile,'Endpoint')== 0 %the left endpoint of the interval is zero
      k=n/2;
    end
    
    x1=linspace(0,1,x1);%computational mesh
    N=length(x1)-1;
    
    
    ordnung=feval(bvpfile,'ordnung'); %vector containign the orders of the solution components
    
    
  elseif (feval(bvpfile,'Infinite') && feval(bvpfile,'EVP'))
    
    %eigenvalue problem on a semi-finite interval
    
    if isempty(werte)
      
      x1=feval(bvpfile,'x1');
      N=length(x1)-1;
      n=feval(bvpfile,'n');
      ordnung=feval(bvpfile,'ordnung');
      
      
    else
      
      
      n=feval(bvpfile,'n');
      
      if feval(bvpfile,'Endpoint') == 0
        k=n/2;
        
        
        %***********************************************************
        %initial profile for the solution component corresponding to
        %the eigenvalue
        %***********************************************************
        werte(end+1,:)=lambda*ones(1,size(werte,2));
        
        %***********************************************************
        %initial profile for the solution component corresponding
        %to the normalization condition (constant function 1)
        %***********************************************************
        werte_end=interp1([0,stellen(end)],[0,1],stellen,'linear','extrap');
        werte(end+1,:)=werte_end;
        
      else
        k=n;
        werte_end=interp1([feval(bvpfile,'Endpoint'),stellen(end)],[0,1],stellen,'linear','extrap');
        werte(end+1,:)=lambda*ones(1,size(werte,2));
        werte(end+1,:)=werte_end;
      end
      
      
      x1=linspace(0,1,x1);
      N=length(x1)-1;
      ordnung=feval(bvpfile,'ordnung');
    end
    
  elseif feval(bvpfile,'EVP') && feval(bvpfile,'Infinite')==0
    
    if isempty(werte)
      
      x1=feval(bvpfile,'x1');
      N=length(x1)-1;
      n=feval(bvpfile,'n');
      
      ordnung=feval(bvpfile,'ordnung');
      
      
    else
      
      
      n=feval(bvpfile,'n');
      ordnung=feval(bvpfile,'ordnung');
      k=n;
      
      
      werte_end=interp1([x1(1),x1(end)],[0,1],stellen,'linear','extrap')  ;
      
      werte(end+1,:)=lambda*ones(1,size(werte,2));
      
      werte(end+1,:)=werte_end;
      
      
    end
    
  elseif feval(bvpfile,'EVP')== 0 && feval(bvpfile,'Infinite')==0
    
    if isempty(werte)
      
      x1=feval(bvpfile,'x1');
      N=length(x1)-1;
      n=feval(bvpfile,'n');
      
      ordnung=feval(bvpfile,'ordnung');
      
      
    else
      
      n=feval(bvpfile,'n'); %number of solution components
      ordnung=feval(bvpfile,'ordnung');
      k=n;
      
    end
    
  end
end

% Raise number of points to minimum specified in options
load options n0; n0=max(n0,8);
if length(x1) < n0
  x1 = pchip(linspace(0,1,length(x1)),x1,linspace(0,1,n0));
  N = length(x1)-1;
end

m=feval(bvpfile,'m');%number of coll.points
parameter=feval(bvpfile,'parameter');

%dummystart=ones(n*N*(m+2),1);
dummystart=ones(N*(n*m+sum(ordnung))+parameter,1);
%equations benoetigt einen Startwert der richtigen Dimension, sonst
%Fehlermeldung
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


if errest_call
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
  
  
else
  
  if ~isempty(werte)
    
    %removal of entries that appear twice in "stellen"
    jj=0;
    for ii=2:length(stellen)
      if stellen(ii)~=stellen(ii-1)
        jj=jj+1;
        stellen2(jj)=stellen(ii-1);
        werte2(1:k,jj)=werte(1:k,jj);
      end
    end
    
    
    stellen2(jj+1)=stellen(end);
    
    werte2(1:k,jj+1)=werte(1:k,end);
  end
  
  
end


for p=1:n
  ordp=ordnung(p);
  interpolationsstellen=zeros(1,(N-1)*(ordp+m-1)+1);
  
  
  if errest_call
    
    interpolationsstellen=zeros(1,(N-1)*(ordnung(p)+m-1)+1);
    for i=0:N-1
      interpolationsstellen(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
    end
    interpolationsstellenwerte=interp1(stellen2,werte2(p,:),interpolationsstellen,'spline');
    
    
    
  else
    
    
    if length(werte) >0
      
      
      
      %the Matlab built in function 'interp1' interpolates and evaluates the
      %underlying function given by werte2
      
      
      
      if feval(bvpfile,'Infinite')
        
        if feval(bvpfile,'Endpoint')== 0
          
          if p>k %components on [1,infty)
            
            
            
            
            
            
            
            ind_1=find(stellen2>=1);
            
            if length(ind_1)<=1
              
              err('initial2');
              return;
              
            end
            
            %transformation to the interval [0,1]
            
            %********************************************************
            %define a computational grid according to the collocation
            %rules
            %********************************************************
            for i=0:N-1
              
              %constructs grid (=mesh + coll.points)
              %inserts coll. points between elements of x1
              interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
              
            end
            
            %****************
            %Interpolation
            %***************
            
            if ind_1(1)-1 <=0
              err('initial1');
              return;
              
            end
            
            if abs(stellen2(ind_1(1)-1)-0)< 10^-4
              stellen2(ind_1(1)-1)=10^-4;
            end
            
            interpolationsstellen=[interpolationsstellen,1./stellen2(ind_1(1)-1)];
            interpolationsstellenwerte=interp1([1./stellen2(ind_1(1)-1),1./stellen2(ind_1)],[werte2(p-k,ind_1(1)-1),werte2(p-k,ind_1)],interpolationsstellen,'linear','extrap');
            
            
            
            
          else
            %if the number of the solution component is <=n/2
            %then no transformation has to be carried out
            
            
            %
            ind_1=find(stellen2<=1);
            
            
            if length(ind_1)<=1
              err('initial1');
              return;
            end
            
            
            for i=0:N-1
              
              
              %construction of the computational grid (=mesh + coll.points)
              %inserts coll. points between elements of x1
              interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
              
            end
            
            
            if ind_1(end)+1 > length(stellen2)
              
              err('initial2');
              return;
              
            end
            
            interpolationsstellen=[interpolationsstellen,stellen2(ind_1(end)+1)];
            interpolationsstellenwerte=interp1([stellen2(ind_1),stellen2(ind_1(end)+1)],[werte2(p,ind_1),werte2(p,ind_1(end)+1)],interpolationsstellen,'spline');
            
            %figure
            %plot(interpolationsstellen,interpolationsstellenwerte,'r')
            
          end
          
          
        else
          
          %the given problem is posed on [a,infty) a ~= 0
          
          
          ind_1=find(stellen2>=1);
          
          
          %transformation to the interval [0,1]
          
          %********************************************************
          %define a computational grid according to the collocation
          %rules
          %********************************************************
          for i=0:N-1
            
            %constructs grid (=mesh + coll.points)
            %inserts coll. points between elements of x1
            interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
            
          end
          
          %****************
          %Interpolation
          %***************
          
          
          
          
          
          interpolationsstellenwerte=interp1(1./stellen2(ind_1),werte2(p,ind_1),interpolationsstellen,'linear','extrap');
          
          
          
          
        end
      else
        %if the problem is posed on a finite interval
        
        for i=0:N-1
          
          % x1 ... new mesh
          %constructs grid (=mesh + coll.points)
          %the number of added collocation points varies with the order
          %of derivative of the solution component z_i
          interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
          
        end
        werte2(p,:);
        interpolationsstellenwerte=interp1(stellen2,werte2(p,:),interpolationsstellen,'spline');
        
      end
      %%%%%%%%%%%%
    else
      
      
      
      
      %************************************************************
      %In case of an eigenvalue problem, the user can insert a guess for
      %the eigenvalue without providing an initial profile for the
      %eigenfunction.
      %In that case the initial profile for the solution components corresponding to the eigenvalue
      %are set to the constant function of the size of the eigenvalue and the
      %others are set to the constant function 1.
      %
      %For a finite interval problem, where no initial profile is provided
      %also the constant function with value 1 serves as starting guess
      %************************************************************
      
      for i=0:N-1
        
        %constructs grid (=mesh + coll.points)
        %the number of added collocation points varies with the order
        %of derivative of the solution component z_i
        interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
        
      end
      
      
      
      if feval(bvpfile,'Infinite')
        if p==n-1 || p==n/2-1
          
          
          if isempty(lambda)
            interpolationsstellenwerte=ones(1,length(interpolationsstellen));
          else
            interpolationsstellenwerte=lambda*ones(1,length(interpolationsstellen));
          end
        else
          
          interpolationsstellenwerte=ones(1,length(interpolationsstellen));
          
        end
        
      else
        
        if p==n-1
          
          if isempty(lambda)
            interpolationsstellenwerte=ones(1,length(interpolationsstellen));
          else
            interpolationsstellenwerte=lambda*ones(1,length(interpolationsstellen));
          end
          
        else
          
          interpolationsstellenwerte=ones(1,length(interpolationsstellen));
        end
        
      end
      
      
      
    end
    
    
  end
  
  
  %***********************************************
  %calculates coefficients in Runge-Kutta basis!!
  %***********************************************
  
  
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
          A(j,sp)= polyval(psi0(sp,:),(intervallstellen(j)-x1((i)+1))/h((i)+1));
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
        j=j+1;
        start(j,1)=intervallstart((i)+1,k,p);
      end
    end
  end
  for k=1:m
    for p=1:n
      j=j+1;
      start(j,1)=intervallstart((i)+1,k+ordnung(p),p);
    end
  end
end
for pii=1:parameter
  j=j+1;start(j,1)=par(pii);
end
start=start(1:j,1);
