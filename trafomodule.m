function [antwort2,antwort6,antwort7] = trafomodule(antwort6,antwort7,anz_glei,antwort_EVP,antwort2,antwort_endpoint,parameter) 
%%*************************************************************************
% TRAFOMODULE transforms a problem posed on [a,inf) (a>=0) to a finite
% interval
% 
% 
% FUNCTION CALL: [antwort2,antwort6,antwort7] = trafomodule(antwort6,antwort7,anz_glei,antwort_EVP,antwort2,antwort_endpoint) 
% INPUT:   antwort6         ... string describing equations (given by user
%                               in GUI)
%          antwort7         ... string describing boundary conditions 
%          anz_glei         ... number of equations
%          antwort_EVP      ... boolean; determines if it is an eigenvalue problem
%          antwort2         ... array with orders of equations
%          antwort_endpoint ... if problem is semi-finite; antwort_endpoint is the left endpoint 
%
%OUTPUT:   antwort2          ... contains orders of new system of equations
%          antwort6         ...  string; describing new system of
%                                equations
%          antwort7         ...  string; describing new boundary conditions
%         
%
% AUTHOR: csimon
% DATE: 02/09
% COMMENT: 
%**************************************************************************


    %******************************************************************** 
    %The transformation is carried out by a change of the independent
    %variable 't'.
    %******************************************************************** 
    
    
    %******************************************************************** 
    %Replacement of all 't' in the original equations
    %by '1/t'. In case that the endpoint 'a' is unequal 0 't' is replaced 
    %by 'a/t'.
    %******************************************************************** 
    ord1=str2num(antwort2);
    
    
    if antwort_endpoint==0
    
       antwort6_brac=strrep(antwort6,']','');

        if strcmp(antwort6_brac(end),';')==0

            antwort6_brac=strcat(antwort6_brac,';');
            antwort6=strcat(antwort6_brac,']');

        end  
        
        equ = regexprep(antwort6,'(?<=^|[^a-zA-Z0-9_])t(?=$|[^a-zA-Z0-9_])','(1/t)');
    else
        equ = regexprep(antwort6,'(?<=^|[^a-zA-Z0-9_])t(?=$|[^a-zA-Z0-9_])',strcat('(',antwort_endpoint,'/t)'));
    end
    
    t = str2sym('t');
    for i=1:anz_glei
        variable_base=strcat('z',num2str(i));
        variable=strcat('z',num2str(i),'(t)');
        
        dashes = '';
        lhs = str2sym(variable);
        for o=1:ord1(i)
            dashes = strcat(dashes,'''');
            equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])',variable_base,dashes,'(?=$|[^''a-zA-Z0-9_])'),strcat(variable_base,'d',num2str(o)));
            lhs = -t^2*diff(lhs,t);
            equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])',variable_base,'d',num2str(o),'(?=$|[^a-zA-Z0-9_])'),strcat('(',strrep(char(lhs),'$','\$'),')'));
        end
        for o=1:ord1(i)
            equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])diff\(z',num2str(i),'\(t\)(,\s*t){',num2str(o),'}\)(?=$|[^a-zA-Z0-9_])'),strcat(variable_base,'d',num2str(o)));
            equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])diff\(z',num2str(i),'\(t\),`\$`\(t,',num2str(o),'\)\)(?=$|[^a-zA-Z0-9_])'),strcat(variable_base,'d',num2str(o)));
        end
        equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])',variable_base,'d2','(?=$|[^''a-zA-Z0-9_])'),strcat(variable_base,''''''));
        equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])',variable_base,'d1','(?=$|[^''a-zA-Z0-9_])'),strcat(variable_base,''''));
    end
    
    %*********************************************************************   
    %Double the number of equations if the problem is posed on [0;infinity).
    %*********************************************************************   
    if antwort_endpoint == 0
     for i=1:anz_glei
        equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])z',num2str(i),'(?=($|d\d+|[^a-zA-Z0-9_]))'),strcat('z',num2str(i+anz_glei)));
     end
     equ=strrep(equ,'[','');
     antwort6=strcat(strrep(antwort6,']',''),equ);
    else
     antwort6 = equ;
    end    
        
    
if strcmp(antwort_EVP,'true')     
  
 ord=str2num(char(antwort2));
 antwort2=strrep(strrep(antwort2,'[',''),']','');
 
 if antwort_endpoint ==0
 
 antwort2=strcat('[',char(antwort2),',',char(antwort2),']');

 
 end 

else
    
 ord=str2num(char(antwort2));
 antwort2=strrep(strrep(antwort2,'[',''),']','');


 if antwort_endpoint ==0
 
 antwort2=strcat('[',char(antwort2),',',char(antwort2),']');
 
 
 end 
 
end
 
 %********************************************************************
 % Construct the Boundary conditions 
 % Note: Only Dirichlet boundary conditions can be posed at infinity!!
 %********************************************************************

 antwort7 = strrep(antwort7,']',''); 

    if strcmp(antwort7(end),';')==0

            antwort7=strcat(antwort7,';');
           
    end

    equ = antwort7;
    if antwort_endpoint==0
        for i=1:anz_glei
            variable_base = strcat('z',num2str(i));
            
            % mapping infinity to 0
            equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])',variable_base,'\(b\)'),strcat('z',num2str(i+anz_glei),'(a)'));

            % continuity
            equ = strcat(equ,variable_base,'(b)=z',num2str(i+anz_glei),'(b);');
            lhs = str2sym(strcat('z',num2str(i+anz_glei),'(t)'));
            for o=1:ord1(i)-1
                lhs = -t^2*diff(lhs,t);
                equ = strcat(equ,variable_base,'d',num2str(o),'(b)=',strcat('(',strrep(char(lhs),'$','\$'),');'));
            end
            for o=1:ord1(i)
                equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])diff\(z',num2str(i+anz_glei),'\(t\)(,\s*t){',num2str(o),'}\)(?=$|[^a-zA-Z0-9_])'),strcat('z',num2str(i+anz_glei),'d',num2str(o),'(b)'));
                equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])diff\(z',num2str(i+anz_glei),'\(t\),\$\(t,',num2str(o),'\)\)(?=$|[^a-zA-Z0-9_])'),strcat('z',num2str(i+anz_glei),'d',num2str(o),'(b)'));
            end
        end
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])t(?=$|[^''a-zA-Z0-9_])','1');
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)d2(?=$|[^''a-zA-Z0-9_])','$1''''');
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)d1(?=$|[^''a-zA-Z0-9_])','$1''');
    else
        % mapping infinity to _tmp
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)\(b\)','$1_tmp(a)');
        
        %mapping a to b
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+d\d+)\(a\)','$1(b)');
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)''''','(?=$|[^''a-zA-Z0-9_])','$1d2');
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)''','(?=$|[^''a-zA-Z0-9_])','$1d1');
        for i=1:anz_glei
            variable_base = strcat('z',num2str(i));
            
            lhs = strcat(variable_base,'(t)');            
            for o=1:ord1(i)
                lhs = -t^2/antwort_endpoint*diff(lhs,t);
                equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])',variable_base,'d',num2str(o),'\(a\)'),strcat('(',strrep(char(lhs),'$','\$'),')'));
            end
            for o=1:ord1(i)
                equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])diff\(z',num2str(i),'\(t\)(,\s*t){',num2str(o),'}\)(?=$|[^a-zA-Z0-9_])'),strcat(variable_base,'d',num2str(o),'(b)'));
                equ = regexprep(equ,strcat('(?<=^|[^a-zA-Z0-9_])diff\(z',num2str(i),'\(t\),\$\(t,',num2str(o),'\)\)(?=$|[^a-zA-Z0-9_])'),strcat(variable_base,'d',num2str(o),'(b)'));
            end
        end
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])t(?=$|[^a-zA-Z0-9_])','1');
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)d2(?=$|[^''a-zA-Z0-9_])','$1''''');
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)d1(?=$|[^''a-zA-Z0-9_])','$1''');

        % mapping _tmp to a
        equ = regexprep(equ,'(?<=^|[^a-zA-Z0-9_])(z\d+)_tmp\(a\)','$1(b)');
    end
    antwort7 = strcat(equ,']');