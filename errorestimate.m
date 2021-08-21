function [x1tau,valerror,maxfehlerwerte,koeff_2,x1tau_2,wert_2x1undtau_2,x1_2,polynome_2,wert_2x1_2,valerror2,tau_infinite,error_infinite,error_infinite2,sol_infinite,lambda,eigenfunction]=errorestimate(coeff,bvpfile,zeichnung,x1,bvpopt,ausgabe,praediktor)

                                                                                                                        %Funktion zum Schätzen des Fehlers einer berechneten Lösung

error_infinite=[];
error_infinite2=[];

if (nargin<7) praediktor=[]; if (nargin<6) ausgabe=1; if (nargin<5) bvpopt=[]; if(nargin<4) x1=[];if(nargin<3) zeichnung=false;end;end;end;end;end;

if (ausgabe==1) || (ausgabe==2)
    fprintf('Initialize error estimate ...\n');
end

%Definiere Werte der zu schätzenden Lösung

if(length(x1)==0)
    x1=feval(bvpfile,'x1');
end
parameter=feval(bvpfile,'parameter');
if parameter>0
    par=coeff(length(coeff)-parameter+1:length(coeff));
else
    par=[];
end
x1tau = equations('x1tau',coeff,bvpfile,x1);
valx1tau = equations('valx1tau',coeff,bvpfile,x1);
standard=feval(bvpfile,'standard');
rho=feval(bvpfile,'rho');
if standard
    m=rho(2);
else
    m=length(rho);
end

% Der Vektor x1 wird verdoppelt (immer an den Mittelpunkten der Intervalle)
% "h/2"

x1_2(1:2:2*length(x1)) = x1;
x1_2(2:2:2*length(x1)-1) = (x1(2:length(x1))+x1(1:length(x1)-1))/2;

%Verwende als Startprofil für die feinere Lösung die zu schätzende Lösung
%Alte Startprofilwerte
[x1_2,start]=initialmesh(bvpfile,x1tau,valx1tau,x1_2,par,[],1);
%Alte Startprofilkoeffizienten
%[x1_2,start]=initialmesh(bvpfile,x1,[],x1_2,[],coeff);


if length(praediktor)>0 %&& length(praediktorlast)==0
    praediktor(3,1)=par(1);
end

%Berechne die feinere Lösung
zeichnung1=0;
[koeff_2,x1_2,wert_2x1_2,x1tau_2,wert_2x1undtau_2,polynome_2,sol_infinite,tau_infinite,lambda,eigenfunction]=run(bvpfile,zeichnung1,x1_2,start,bvpopt,ausgabe,praediktor);


% Werte die feinere Lösung an den Stellen der gröberen aus

wert_2x1undtau = equations('wert',koeff_2,bvpfile,x1_2,x1tau);

valerror = (2^m / (1-2^m)) * (wert_2x1undtau - valx1tau);


% Werte die gröbere Lösung an den Stellen der feineren aus (für die
% Umkehrung)

wert_1x1undtau_2 = equations('wert',coeff,bvpfile,x1,x1tau_2);

valerror2= (1 / (1-2^m)) * (wert_2x1undtau_2 - wert_1x1undtau_2);



testfehlerwerte2=wert_2x1undtau_2 - wert_1x1undtau_2;


% Definition der Maximumnorm der Absolutbeträge der Fehlerwerte

if length(valerror(:,1))>1
    maxfehlerwerte=max(abs(valerror));
else
    maxfehlerwerte=abs(valerror);
end

if zeichnung
    
    figure;
    n=feval(bvpfile,'n');
    


%bvpfile=strcat(pfad,bvpfile)



    
    
    
    if feval(bvpfile,'Infinite') %posed on infinite interval
        
        if feval(bvpfile,'EVP')  %eigenvalue problem 
         
            %valerror
           if feval(bvpfile,'Endpoint')==0  
               
              
               for i=1:length(valerror(:,1))/2-2
                 
                   subplot(length(valerror(:,1))/2-2,1,i);

                   [error_infinite(i,:),tau_infinite]=backtransf(valerror(i,2:end),valerror(i+n/2,2:end),x1tau(2:end),0);
                  
                   plot(tau_infinite,real(error_infinite(i,:)),'color','black');
                   
                   xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
               end
            
                    title( subplot(length(valerror(:,1))/2-2,1,1),'Absolute error on the coarse grid');
               
               
           else %endpoint is not 0

                for i=1:length(valerror(:,1))-2
                  
                  subplot(length(valerror(:,1))-2,1,i);

                  [error_infinite(i,:),tau_infinite]=backtransf([],valerror(i,2:end),x1tau(2:end),feval(bvpfile,'Endpoint'));
                 
                   title('Absolute error on the coarse grid');
                xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                end 
                
                plot(subplot(length(valerror(:,1))-2,1,1),tau_infinite,real(error_infinite(i,:)),'color','black'); 
           end 
           %valerror2
            
            figure;
           
            if feval(bvpfile,'Endpoint')==0
                for i=1:length(valerror2(:,1))/2-2
                    
                    subplot(length(valerror2(:,1))/2-2,1,i);
                    [error_infinite2(i,:),tau_infinite]=backtransf(valerror2(i,2:end),valerror2(i+n/2,2:end),x1tau_2(2:end),0);
                   plot(tau_infinite,real(error_infinite2(i,:)),'color','black')  
                   
                  xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);

                end    
                
                title(subplot(length(valerror2(:,1))/2-2,1,1),'Absolute error on the coarse grid of the solution from the fine grid','FontSize',9);


            else
                for i=1:length(valerror2(:,1))-2
                  
                    subplot(length(valerror2(:,1))-2,1,i);
                  [error_infinite2(i,:),tau_infinite]=backtransf([],valerror2(i,2:end),x1tau_2(2:end),feval(bvpfile,'Endpoint'))  ;
                  plot(tau_infinite,real(error_infinite2(i,:)),'color','black');
                 
                xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                end 
                
                 title(subplot(length(valerror2(:,1))-2,1,1),'Absolute error on the coarse grid of the solution from the fine grid','FontSize',9);
                
                
            end 
            
        
        else % no eigenvalue problem

            %valerror
            if feval(bvpfile,'Endpoint')==0 
               for i=1:length(valerror(:,1))/2
                   
                   subplot(length(valerror(:,1))/2,1,i);

                   [error_infinite(i,:),tau_infinite]=backtransf(valerror(i,2:end),valerror(i+n/2,2:end),x1tau(2:end),0);
                   plot(tau_infinite,real(error_infinite(i,:)),'color','black')
                     title('Absolute error on the coarse grid');
                   xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                 
               end
            
               
               
            else %endpoint is not 0
            
                for i=1:length(valerror(:,1))
                  
                  subplot(length(valerror(:,1)),1,i);
                  [error_infinite(i,:),tau_infinite]=backtransf([],valerror(i,2:end),x1tau(2:end),feval(bvpfile,'Endpoint'));
                  plot(tau_infinite,real(error_infinite(i,:)),'color','black'); 
                  
                  xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                end 
                  title( subplot(length(valerror(:,1)),1,1),'Absolute error on the coarse grid');
               
           end 
            
            figure;
           %valerror2
            if feval(bvpfile,'Endpoint')==0
                for i=1:length(valerror2(:,1))/2
                    
                 
                   [error_infinite2(i,:),tau_infinite]=backtransf(valerror2(i,2:end),valerror2(i+n/2,2:end),x1tau_2(2:end),0);
                   plot(tau_infinite,real(error_infinite2(i,:)),'color','black');         
                  
                   xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                end    
                
                 title(subplot(length(valerror2(:,1))/2,1,1),'Absolute error on the coarse grid of the solution from the fine grid','FontSize',9);

            else %endpoint is not 0
                for i=1:length(valerror2(:,1))
                 
                  subplot(length(valerror2(:,1)),1,i);
                  [error_infinite2(i,:),tau_infinite]=backtransf([],valerror2(i,2:end),x1tau_2(2:end),feval(bvpfile,'Endpoint'))  ;
                  plot(tau_infinite,real(error_infinite2(i,:)),'color','black');
                 
                  xlabel('t');
                  ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                end 
                  title( subplot(length(valerror2(:,1)),1,1),'Absolute error on the coarse grid of the solution from the fine grid','FontSize',9);
            end 
            
        


        end
        
    else %problem is posed on a finite interval 

       
        
        if feval(bvpfile,'EVP')
            
            
            for i=1:length(valerror(:,1))-2
                
                subplot(length(valerror(:,1))-2,1,i);
                plot(x1tau,real(valerror(i,:)),'color','black');
                xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
            end
            
             title(subplot(length(valerror(:,1))-2,1,1),'Absolute error on the coarse grid');
            
            figure;
             
            for i=1:length(valerror2(:,1))-2
                
                subplot(length(valerror2(:,1))-2,1,i);
                plot(x1tau_2,real(valerror2(i,:)),'color','black');
                xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
            end           
            
            
            title(subplot(length(valerror2(:,1))-2,1,1),'Absolute error on the coarse grid of the solution from the fine grid','FontSize',9);
            
        else %no eigenvalue problem
             
            for i=1:length(valerror(:,1))
                
                subplot(length(valerror(:,1)),1,i);
                plot(x1tau,real(valerror(i,:)),'color','black');
                xlabel('t');
                ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                
            end
            
            title(subplot(length(valerror(:,1)),1,1),'Absolute error on the coarse grid');
            figure;
             
            for i=1:length(valerror2(:,1))
                
                subplot(length(valerror2(:,1)),1,i);
                plot(x1tau_2,real(valerror2(i,:)),'color','black');
                
                 xlabel('t');
                 ylabel(['$\Delta_{z_',num2str(i),'}$'],'interpreter','latex','FontSize',14);
                
            end
            
            title(subplot(length(valerror2(:,1)),1,1),'Absolute error on the coarse grid of the solution from the fine grid','FontSize',9);
        end 
    end 
    
end


% This one is used only here - later initialmesh should serve all purposes
% function [x1,start]=initialmesh2(bvpfile,stellen,werte,x1,par)
% if (nargin<6) koeff=[]; if (nargin<5) par=[]; end;end;
% 
% N=length(x1)-1;
% n=feval(bvpfile,'n');
% m=feval(bvpfile,'m');
% ordnung=feval(bvpfile,'ordnung');
% parameter=feval(bvpfile,'parameter');
% dummystart=ones(N*(n*m+sum(ordnung))+parameter,1);
% faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];
% if length(par)==0
%     par=ones(parameter);
% end
% if length(koeff)>0 && parameter>0
%     par=koeff(end-parameter+1:end);
% end
% if max(ordnung)>0
%     psi=equations('basispolynome',dummystart,bvpfile,x1);
% end
% if min(ordnung)==0
%     psi0=equations('basispolynome0',dummystart,bvpfile,x1);
% end
% for i=0:N-1
%     h(i+1)=x1(i+2)-x1(i+1);
% end
% 
% 
% for ord=1:max(ordnung)
%     proj_intervallstellen=0:1/(m+ord-1):1;
%     for i=1:m
%         psival(ord,i,1:m+ord)=polyval(psi(i,:,ord),proj_intervallstellen);
%     end
% end
% jj=0;
% for ii=2:length(stellen)
%     if stellen(ii)~=stellen(ii-1)
%         jj=jj+1;
%         stellen2(jj)=stellen(ii-1);
%         werte2(1:n,jj)=werte(1:n,jj);
%     end
% end
% stellen2(jj+1)=stellen(end);
% werte2(1:n,jj+1)=werte(1:n,end);
% for p=1:n
%     ordp=ordnung(p);
%     interpolationsstellen=zeros(1,(N-1)*(ordnung(p)+m-1)+1);
%     for i=0:N-1
%         interpolationsstellen(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
%     end
%     interpolationsstellenwerte=interp1(stellen2,werte2(p,:),interpolationsstellen,'spline');
%   
%     %until here
%     for i=0:N-1
%         if (m==1 && ordnung(p)==0)
%             intervallstellen(1)=x1((i+1));
%         else
%             intervallstellen=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
%         end
%       % if length(koeff)==0
%             %intervallwerte=interp1(stellen,werte(p,:),intervallstellen,'spline');
%             intervallwerte=interpolationsstellenwerte(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1);
%       %  else
%       %      intervallwerte=equations('wert',koeff,bvpfile,stellen,intervallstellen);
%       %      intervallwerte=intervallwerte(p,:);
%       %  end
%         
%         A=zeros(m+ordnung(p),ordnung(p)+m);
%         for j=1:m+ordnung(p)
%             for q=1:ordnung(p)
%                 A(j,q)=(intervallstellen(j)-x1((i)+1))^(q-1)/faktorielle(q);
%             end
%             
%             
%             for sp=1:m
%                 
%                 if ordp>0
%                     %A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*polyval(psi(sp,:,ordnung(p)),(intervallstellen(j)-x1(i+1))/h(i+1));
%                     A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*psival(ordnung(p),sp,j);
%                 else                    
%                     A(j,sp)=polyval(psi0(sp,:),(intervallstellen(j)-x1((i)+1))/h((i)+1));
%                 end
%             end
%         end
%         help=A\intervallwerte';
%         intervallstart((i)+1,1:length(help),p)=help;
%     end
% end
% 
% j=0;
% start=zeros(N*max(ordnung)*n+N*m*n+parameter,1);
% for i=0:N-1
%     for k=1:max(ordnung)
%         for p=1:n
%             if k<=ordnung(p)
%                 j=j+1;start(j,1)=intervallstart((i)+1,k,p);
%             end
%         end
%     end
%     for k=1:m        
%         for p=1:n
%            j=j+1;start(j,1)=intervallstart((i)+1,k+ordnung(p),p);
%         end
%     end
% end
% for pii=1:parameter
%     j=j+1;start(j,1)=par(pii);
% end
% start=start(1:j,1);
