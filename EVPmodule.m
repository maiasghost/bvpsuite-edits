function [antwort2,antwort6,antwort7] = EVPmodule(antwort2,antwort6x,antwort7x)  
%************************************************************************
% EVPMODULE recasts an Eigenvalue Problem into a Boundary Value Problem
% by adding two additional equations and boundary conditions
% 
% 
%FUNCTION CALL:
%[antwort2,antwort6,antwort7] = EVPmodule(antwort2,antwort6x,antwort7x)  
%INPUT:    antwort6x         ...  string,  contains the problem defining equations (given by user
%                                in GUI);
%          antwort7x         ...  string, which contains the boundary
%                                conditions of the given problem;
%          antwort2          ... array, which stores the orders of equations
%
%OUTPUT:   antwort2          ... contains orders of 'augmented' system of equations
%          antwort6         ... string; constitutes 'augmented' system of
%                                equations
%          antwort7         ... string; contains 'augmented' boundary conditions

%         
%
% AUTHOR: csimon
% DATE: 02/09
% COMMENT: for the mathematical background see "Collocation Methods for the
% Solution of Eigenvalue Problems for Singular Ordinary Differential
% Equations" by Auzinger, Karner, Koch and Weinm|ller 
%**************************************************************************



  ord=str2num(char(antwort2));  
  %Since 2 solution components are added, the corresponding orders of
  %derivatives have to be inserted.
  ord_new = strcat('[',char(antwort2),',1,1]'); 
  antwort2= ord_new;
  ord_new= str2num(ord_new);
 
  
   antwort6 = strrep(antwort6x,']','');
   


   hilfsvar=num2str(length(ord)+1);  
   hilfsvar=strcat('z',hilfsvar);
   
   %**********************************************
   %**********************************************
   %Adding the new solution components:
   %**********************************************
   %**********************************************

   
   %Replacement of the string 'lambda' by
   %z_n+1; n=length(ord);
   antwort6=strrep(antwort6,'lambda',hilfsvar);
   
   %*****************************************************
   %Insertion of the new additional equation z_{n+1}'=0 :
   %*****************************************************
   hilfsvar=strcat(hilfsvar,'''=0;');
   antwort6=strcat(char(antwort6),';',hilfsvar);
     
   hilfsvar=strcat('z',num2str(length(ord)+2),'''');
   
   
   %*****************************************************
   %Insertion of the new additional equation
   % z_{n+2}'= z_1^2+z_2^2+...+z_n^2:
   %*****************************************************
   
   hilfsvar2=strcat('z',num2str(1),'^2');
   
   
   antwort6_new = regexp(antwort6,'[^;]*','match');
   
   
     
      c=antwort6_new{1};
     
      
   
        for i=2:length(ord)
   
        
            hilfsvar2 = strcat(hilfsvar2,'+','z',num2str(i),'^2');
     
        end 
   
     
   
   
   antwort6=strcat(antwort6,hilfsvar,'=',hilfsvar2,';');
   antwort6= strcat(char(antwort6),']');

   
   %*****************************************************
   %        Construction of the Boundary Conditions 
   %*****************************************************
   %Insertion of the new boundary conditions associated with 
   %the additonal equations.
   %There are no boundary conditions added for the trivial equation
   %lambda'=0;
   
   antwort7_old=antwort7x;
   antwort7 = strrep(antwort7x,']',''); 
   
   
        if strcmp(antwort7(end),';')==0

            antwort7=strcat(antwort7,';');
           
        end


   antwort7=strcat(antwort7,strcat('z',num2str(length(ord)+2),'(a)','=0;'),strcat('z',num2str(length(ord)+2),'(b)','=1;'));
   antwort7= strcat(char(antwort7),']');
