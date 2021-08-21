function plot_results(t,values,n,EVP,Infinite,left,right)


figure

for i=1:n
    

        
        subplot(n,1,i);
        plot(t,values(i,:),'color','black');
        xlabel('t');
        ylabel(['z_',num2str(i)]);
        
        
     
        
        if Infinite == 1
            xlim([left right])
        end
          
end
   
           
        if EVP==0
            title(subplot(n,1,1),'Numerical solution');
        else
             a=values(n+1,1);
            title(subplot(n,1,1),['Eigenfunction for eigenvalue \lambda = ',num2str(a)]);
        end



   
   