function varargout =evpsolve(varargin)
% EVPSOLVE M-file for evpsolve.fig
%      EVPSOLVE, by itself, creates a new EVPSOLVE or raises the existing
%      singleton*.
%
%      H = EVPSOLVE returns the handle to a new EVPSOLVE or the handle to
%      the existing singleton*.
%
%      EVPSOLVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVPSOLVE.M with the given input arguments.
%
%      EVPSOLVE('Property','Value',...) creates a new EVPSOLVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before evp_fd_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to evpsolve_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help evpsolve

% Last Modified by GUIDE v2.5 21-Sep-2009 14:38:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @evpsolve_OpeningFcn, ...
    'gui_OutputFcn',  @evpsolve_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
%if nargin && ischar(varargin{1})
%   gui_State.gui_Callback = str2func(varargin{1});
%end

if nargin && ischar(varargin{1})


    if strfind(varargin{1},'.m')==1
        filename=varargin{1};
    else


        gui_State.gui_Callback = str2func(varargin{1});
    end
end



if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%set([handles.calculate_button,handles.eigenvector_button,handles.plot_button,...
%        handles.close_button,handles.save_button,handles.help_button,handles.selectAll_button],'Enable','off');


% --- Executes just before evpsolve is made visible.
function evpsolve_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to evpsolve (see VARARGIN)
% Choose default command line output for evpsolve
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes evpsolve wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.printev,'Value',1);



% --- Outputs from this function are returned to the command line.
function varargout = evpsolve_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Specify the function ''g'' with the variables z_1(t),...,z_n(t),z_1''(t),...,z_n''(t).The equations have to be posed in the form g(t,z1,...,zn,z1'',...,zn'') = lambda*f(t,z1,...,zn). Denote the eigenvalue with lambda. Note that only first order equations can be solved. Please type always ''1*z_i'' if z_i appears without coefficient in the equation. ','Help');


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Specify the additional/boundary conditions seperated by semicolon. (Number: Sum of orders and parameters). z_i(a) =: zi(a), z_i''(a) =:zi''(a), z_i(b) =: zi(b), z_i''(b) =: zi''(b). z and z1 are equal. E.g.: z1(a) = 3;  z2(a) = 0;  z1(b) = 0 / For inner-interval conditions write z2(c1)=...;z1(c3)=... Attention: a and b must not be used!','Help');



function editOrders_Callback(hObject, eventdata, handles)
% hObject    handle to editOrders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOrders as text
%        str2double(get(hObject,'String')) returns contents of editOrders as a double


% --- Executes during object creation, after setting all properties.
function editOrders_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOrders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Specify the orders of the differential equations. The code can only deal with systems of order one or order two, e.g. a system of 2 equations of order two should be denoted by the row vector [2,2]. ','Help');



%%%%%%%%%%%%%%%%%  SAVE     %%%%%%%%%%%%%%
% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read GUI input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mainGUIhandle=bvpsuite;
mainGUIdata  = guidata(mainGUIhandle);


filename=get(mainGUIdata.edit1,'String');

check_filename=strfind(filename,'.m');
check_bvpsuite=strfind(filename,'bvpsuite.m');



if isempty(check_filename) || (isempty(check_bvpsuite)==0)
    err('bvps_errdlg06');
    return; 
end     
pfad=get(mainGUIdata.pfad,'String');



%filename = strcat('evp_',filename);


low_ind=get(handles.low_ind,'String');
up_ind=get(handles.up_ind,'String');

low_ind=str2num(low_ind);
up_ind=str2num(up_ind);

antwort_G =get(handles.editG,'String');
antwort_A =get(handles.editA0,'String');
antwort_Ba =get(handles.editB0,'String');
antwort_Bb =get(handles.editB1,'String');
antwort_A1 =get(handles.editA1,'String');
antwort_orders= get(handles.editOrders,'String'); %orders of ODES

if strcmp(antwort_G,'[]')
   
 errordlg('G can''t be empty. The code can only deal with eigenvalue problems. If you want to solve a BVP use the main GUI. Check your input! ','error');   

 return;   
 
end



%antwort_discret2 = regexp(antwort_discret,'[^,]*','match');

%antwort_infinite is true, if the problem is posed
%on a semi-finite interval
if get(handles.radiobuttonInfinite,'Value')==1
    antwort_infinite='true';

else
    antwort_infinite='false';
end

orders=str2num(char(antwort_orders));



%change boundary conditions if the problem is posed
%on a semi-infinite interval


%antwort_endpoint is empty if no semi-infinite interval problem 
%otherwise it is 0

if strcmp(antwort_infinite,'false')
    a = str2num(get(handles.edita,'String'));
    b = str2num(get(handles.editb,'String'));
    N = str2num(get(handles.editN,'String'));
    antwort_endpoint=[];
else
    
    a = str2num(get(handles.edita,'String'));
    N = str2num(get(handles.editN,'String'));
    antwort_endpoint=a;
    
    if a==0  
        b=1;
    else
        a=0;
        b=1;
    end
    
   
end



%%%%%%%%%%%%%%%%%%%%%



%If the user prefers to type in the matrices to
%specify the equations
%then the matrices are immediately written in the evpfile
n=length(str2num(antwort_orders));
antwort_orders_old=antwort_orders;





antwort_A_old=antwort_A;
antwort_A1_old=antwort_A1;
antwort_G_old=antwort_G;
antwort_Ba_old=antwort_Ba;
antwort_Bb_old=antwort_Bb;

if strcmp(antwort_infinite,'false') 
  
    
    if max(str2num(antwort_orders))==1 %first order problem
    
        
        %if you have a first order problem and a finite interval
        
        
        antwort_Ba=antwort_Ba_old;
        antwort_Bb=antwort_Bb_old;
        antwort_A_old=antwort_A;
        antwort_A1_old=antwort_A1;
        antwort_G_old=antwort_G;
        
        
        
    else    
        
    %if your have a second order problem
    %and finite interval, you have to approximate
    %the boundary conditions involving derivatives
    
    %change boundary conditions:
    antwort_Ba=strrep(antwort_Ba,'[',';');
    antwort_Ba=strrep(antwort_Ba,']',';');
    Ba_rows = regexp(antwort_Ba,'[^;]*','match'); 
    for i=1:length(Ba_rows)

        Ba_col(i,:)= regexp(Ba_rows(i),'[^,]*','match');

    end

    antwort_Bb=strrep(antwort_Bb,'[',';');
    antwort_Bb=strrep(antwort_Bb,']',';');
    Bb_rows = regexp(antwort_Bb,'[^;]*','match');
    for i=1:length(Bb_rows)

        Bb_col(i,:)= regexp(Bb_rows(i),'[^,]*','match');

    end
    
    
    if length(Ba_rows) ~= length(Bb_rows) || size(Bb_col,2) ~= size(Ba_col,2)
        
        errordlg('Check your input! Your boundary conditions might have the wrong size! Make sure that your problem is well-posed.','error');
        return;
    end 

    nbc=length(Ba_rows);

    Ba=cell(nbc,3*n); %you need 3*n rows, because you need y_0,y_1 and y_2 for the approx
                      %of the derivative in the bc
    Bb=cell(nbc,3*n);


    %bcs at a
    for i=1:nbc
  
      
        for j=1:n
            
            Ba{i,j} = strcat(Ba_col{i}(j));
            Ba{i,j}=strcat(Ba{i,j},'+',Ba_col{i}(n+j),'*(-3/2)');
            Ba{i,j+n}=strcat('+',Ba_col{i}(n+j),'*(2)');
            Ba{i,j+2*n}=strcat('+',Ba_col{i}(n+j),'*(-1/2)');


        end


    end

    
    
    %bcs at b
    for i=1:nbc


        for j=1:n

           
            Bb{i,2*n+j} = strcat(Bb_col{i}(j));

            Bb{i,j}=strcat('+',Bb_col{i}(n+j),'*(1/2)');
            Bb{i,j+n}=strcat('+',Bb_col{i}(j+n),'*(-2)');
            Bb{i,j+2*n}=strcat(Bb{i,2*n+j},'+',Bb_col{i}(j+n),'*(3/2)');


        end


    end


    
    
    %create strings to that it can be written into evpfile
    antwort_Ba='[';

    for i=1:size(Ba,1)


        for j=1:size(Ba,2)-1


            if isempty(Ba{i,j})
                antwort_Ba=strcat(antwort_Ba,'0',',');
            else
                antwort_Ba=strcat(antwort_Ba,Ba{i,j},',');
            end

        end

        if isempty(Ba{i,size(Ba,2)})
            antwort_Ba=strcat(antwort_Ba,'0',';');
        else
            antwort_Ba=strcat(antwort_Ba,Ba{i,size(Ba,2)},';');
        end

    end

    antwort_Ba=strcat(antwort_Ba,']');
    antwort_Ba=antwort_Ba{1};


    antwort_Bb='[';

    for i=1:size(Bb,1)


        for j=1:size(Bb,2)-1

            if isempty(Bb{i,j})

                antwort_Bb=strcat(antwort_Bb,'0',',');
            else
                antwort_Bb=strcat(antwort_Bb,Bb{i,j},',');
            end
        end

        if isempty(Bb{i,size(Bb,2)})

            antwort_Bb=strcat(antwort_Bb,'0',';');
        else

            antwort_Bb=strcat(antwort_Bb,Bb{i,size(Bb,2)},';');
        end
    end

    antwort_Bb=strcat(antwort_Bb,']');
    antwort_Bb=antwort_Bb{1};

    end 


else
%**********************************************    
%if your problem is posed on a semi-finite interval
%**********************************************

    %set up the matrices A,G
    antwort_A=strrep(antwort_A,'[',';');
    antwort_A=strrep(antwort_A,']',';');
    A_rows=regexp(antwort_A,'[^;]*','match');


    for i=1:length(A_rows)

        A_col(i,:)= regexp(A_rows(i),'[^,]*','match');

    end

    
    antwort_G=strrep(antwort_G,'[',';');
    antwort_G=strrep(antwort_G,']',';');
    G_rows=regexp(antwort_G,'[^;]*','match');


    for i=1:length(G_rows)

        G_col(i,:)= regexp(G_rows(i),'[^,]*','match');

    end




if antwort_endpoint ==0
   
    %double the number of rows and columns because 1 equation becomes two
    %equ.
    antwort_A=strcat('[');

    for i=1:length(A_rows)



        antwort_A=strcat(antwort_A,A_rows{i});


        for j=1:length(A_col{i})-1


            antwort_A= strcat(antwort_A,',','0');

        end

        antwort_A=strcat(antwort_A,',','0',';');

    end

else
    
    
     antwort_A=strcat('[');

    
    
end 



    for i=1:length(A_rows)


     if antwort_endpoint==0   
        for j=1:length(A_col{i})


            antwort_A= strcat(antwort_A,'0',',');

        end

     end     

        for j=1:length(A_col{i})-1

            if max(str2num(antwort_orders))==2
                
                if antwort_endpoint==0
                    antwort_A=strcat(antwort_A,'(1/t^4)*','(',strrep(A_col{i}{j},'t','(1/t)'),')',',');
                    antwort_A =  strrep(antwort_A,'sqr(1/t)','sqrt');
                    
                else
                    antwort_A=strcat(antwort_A,strcat('(',num2str(antwort_endpoint),'^2/t^4)*'),'(',strrep(A_col{i}{length(A_col{i})},'t',strcat('(',num2str(antwort_endpoint),'/t)')),')',';');
                     antwort_A =  strrep(antwort_A,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
                end
                
            else
                
                if antwort_endpoint==0
                    
                    antwort_A=strcat(antwort_A,'(-1/t^2)*','(',strrep(A_col{i}{j},'t','(1/t)'),')',',');
                    antwort_A =  strrep(antwort_A,'sqr(1/t)','sqrt');
                    
                else
                    antwort_A=strcat(antwort_A,'(-',num2str(antwort_endpoint),'/t^2)*','(',strrep(A_col{i}{j},'t',strcat('(',num2str(antwort_endpoint),'/t)')),')',',');
                    antwort_A =  strrep(antwort_A,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
                    
                end
                
            end
        end

        if max(str2num(antwort_orders))==2    
            if antwort_endpoint==0   
                antwort_A=strcat(antwort_A,'(1/t^4)*','(',strrep(A_col{i}{length(A_col{i})},'t','(1/t)'),')',';');
                antwort_A =  strrep(antwort_A,'sqr(1/t)','sqrt');
    
            else
                antwort_A=strcat(antwort_A,strcat('(',num2str(antwort_endpoint),'^2/t^4)*'),'(',strrep(A_col{i}{length(A_col{i})},'t',strcat('(',num2str(antwort_endpoint),'/t)')),')',';');
                 antwort_A =  strrep(antwort_A,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
            end
        else

            if antwort_endpoint==0 
                antwort_A=strcat(antwort_A,'(-1/t^2)*','(',strrep(A_col{i}{length(A_col{i})},'t','(1/t)'),')',';');
                antwort_A =  strrep(antwort_A,'sqr(1/t)','sqrt');
    
            else
                antwort_A=strcat(antwort_A,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(A_col{i}{length(A_col{i})},'t',strcat('(',num2str(antwort_endpoint),'/t)')),')',';');
                antwort_A =  strrep(antwort_A,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
            
            end
        end    

    end
    antwort_A=strcat(antwort_A,']');



 
    % only if you have a second order problem, you also have to transform the
    %matrix A1
    if max(str2num(antwort_orders))==2


    antwort_A1=strrep(antwort_A1,'[',';');
    antwort_A1=strrep(antwort_A1,']',';');
    A1_rows = regexp(antwort_A1,'[^;]*','match');


    for i=1:length(A1_rows)

        A1_col(i,:)= regexp(A1_rows(i),'[^,]*','match');

    end


    if antwort_endpoint ==0
    antwort_A1=strcat('[');
    for i=1:length(A1_rows)



        antwort_A1=strcat(antwort_A1,A1_rows{i});


        for j=1:length(A1_col{i})-1


            antwort_A1= strcat(antwort_A1,',','0');

        end

        antwort_A1=strcat(antwort_A1,',','0',';');

    end
    else
      antwort_A1=strcat('[');  
        
    end 
   




    for i=1:length(A1_rows)

        if antwort_endpoint ==0

            for j=1:length(A1_col{i})


                antwort_A1= strcat(antwort_A1,'0',',');

            end

        end
        
        for j=1:length(A1_col{i})-1

            if i==j
                
                if antwort_endpoint==0
                antwort_A1=strcat(antwort_A1,'(-1/t^2)*','(',strrep(A1_col{i}{j},'t','(1/t)'),')','-2/t',',');
                 antwort_A1 =  strrep(antwort_A1,'sqr(1/t)','sqrt');
                
                else
                  antwort_A1=strcat(antwort_A1,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(A1_col{i}{j},'t','(1/t)'),')','-2/','(',num2str(antwort_endpoint),'^2*t)',',');   
                  antwort_A1 =  strrep(antwort_A1,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');  
                
                end 
                
                
            else
                
                if antwort_endpoint==0
                antwort_A1=strcat(antwort_A1,'(-1/t^2)*','(',strrep(A1_col{i}{j},'t','(1/t)'),')',',');
                 antwort_A1 =  strrep(antwort_A1,'sqr(1/t)','sqrt');
                else
                 antwort_A1=strcat(antwort_A1,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(A1_col{i}{j},'t','(1/t)'),')',',');   
                  antwort_A1 =  strrep(antwort_A1,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');  
                end 
                
            end

        end

        if i==length(A1_col{i})
            
            if antwort_endpoint==0
            antwort_A1=strcat(antwort_A1,'(-1/t^2)*','(',strrep(A1_col{i}{length(A1_col{i}) },'t','(1/t)'),')','-2/t',';');
              antwort_A1 =  strrep(antwort_A1,'sqr(1/t)','sqrt');
            else
             antwort_A1=strcat(antwort_A1,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(A1_col{i}{length(A1_col{i}) },'t','(1/t)'),')','-2/','(',num2str(antwort_endpoint),'^2*t)',';');
              antwort_A1 =  strrep(antwort_A1,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
            end 
            
        else
            
            if antwort_endpoint==0
            antwort_A1=strcat(antwort_A1,'(-1/t^2)*','(',strrep(A1_col{i}{length(A1_col{i}) },'t','(1/t)'),')',';');
              antwort_A1 =  strrep(antwort_A1,'sqr(1/t)','sqrt');
            else
              antwort_A1=strcat(antwort_A1,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(A1_col{i}{length(A1_col{i}) },'t','(1/t)'),')',';');   
               antwort_A1 =  strrep(antwort_A1,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
            end 
            
        end



    end
    antwort_A1=strcat(antwort_A1,']');
 
   
    end 
    

    if antwort_endpoint ==0
    antwort_G=strcat('[');

    for i=1:length(G_rows)



        antwort_G=strcat(antwort_G,G_rows{i});


        for j=1:length(G_col{i})-1


            antwort_G= strcat(antwort_G,',','0');

        end

        antwort_G=strcat(antwort_G,',','0',';');

    end

    else
      antwort_G=strcat('[');   
        
    end 
    
    

    for i=1:length(G_rows)

     
       if antwort_endpoint ==0 
        for j=1:length(G_col{i})


            antwort_G= strcat(antwort_G,'0',',');

        end
       end 

        for j=1:length(G_col{i})-1

            if max(str2num(antwort_orders))==2
            
                if antwort_endpoint==0
                
                    antwort_G=strcat(antwort_G,'(1/t^4)*','(',strrep(G_col{i}{j},'t','(1/t)'),')',',');
                    antwort_G =  strrep(antwort_G,'sqr(1/t)','sqrt');
                
                else
                    
                antwort_G=strcat(antwort_G,strcat('(',num2str(antwort_endpoint),'^2/t^4)*'),'(',strrep(G_col{i}{j},'t','(1/t)'),')',',');   
                 antwort_G =  strrep(antwort_G,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
                end 
            else
                
            if antwort_endpoint==0    
                antwort_G=strcat(antwort_G,'(-1/t^2)*','(',strrep(G_col{i}{j},'t','(1/t)'),')',',');    
                antwort_G =  strrep(antwort_G,'sqr(1/t)','sqrt');
            else

               antwort_G=strcat(antwort_G,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(G_col{i}{j},'t','(1/t)'),')',',');  
               antwort_G =  strrep(antwort_G,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
            end 
            
            end
            
        end

        if max(str2num(antwort_orders))==2
            
        if antwort_endpoint==0    
        antwort_G=strcat(antwort_G,'(1/t^4)*','(',strrep(G_col{i}{length(G_col{i})},'t','(1/t)'),')',';');
        antwort_G =  strrep(antwort_G,'sqr(1/t)','sqrt');
        else
        antwort_G=strcat(antwort_G,strcat('(',num2str(antwort_endpoint),'^2/t^4)*'),'(',strrep(G_col{i}{length(G_col{i})},'t','(1/t)'),')',';');
           antwort_G =  strrep(antwort_G,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
        end
        else
            
          if antwort_endpoint==0  
             antwort_G=strcat(antwort_G,'(-1/t^2)*','(',strrep(G_col{i}{length(G_col{i})},'t','(1/t)'),')',';');  
             antwort_G =  strrep(antwort_G,'sqr(1/t)','sqrt');
          else
             antwort_G=strcat(antwort_G,strcat('(-',num2str(antwort_endpoint),'/t^2)*'),'(',strrep(G_col{i}{length(G_col{i})},'t','(1/t)'),')',';');    
             antwort_G =  strrep(antwort_G,strcat('sqr(',num2str(antwort_endpoint),'/t)'),'sqrt');
          end 
          
        end


    end
    antwort_G=strcat(antwort_G,']');


    

  

    antwort_Ba=strrep(antwort_Ba,'[',';');
    antwort_Ba=strrep(antwort_Ba,']',';');
    antwort_Ba
    Ba_rows = regexp(antwort_Ba,'[^;]*','match');
    
    
    for i=1:length(Ba_rows)

        Ba_col(i,:)= regexp(Ba_rows(i),'[^,]*','match')
        
       
       if i>1 
           
           if size(Ba_col{i},2) ~= size(Ba_col{i-1},2)
               
                err_evp('evp2');
                return;
               
           end
       end 
    end

    
   size(Ba_col(1,:),1)
   Ba_col(1,:)
       
   size(Ba_col{2},2)
   size(Ba_col{1},2)

    antwort_Bb=strrep(antwort_Bb,'[',';');
    antwort_Bb=strrep(antwort_Bb,']',';');
    Bb_rows = regexp(antwort_Bb,'[^;]*','match');
    for i=1:length(Bb_rows)

        Bb_col(i,:)= regexp(Bb_rows(i),'[^,]*','match');
        
       if i>1 
           
           if size(Bb_col{i},2) ~= size(Bb_col{i-1},2)
               
                err_evp('evp2');
                return;
               
           end
       end 

    end

 
  
    
  if length(Ba_rows) ~= length(Bb_rows) || size(Bb_col,2) ~= size(Ba_col,2) 
        
        err_evp('evp3');
        return;
    end 
    
    
    nbc=length(Ba_rows);
    
if antwort_endpoint==0  

     %if second order problem you have derivatives in bcs 
    
    %nbc defines the number of boundary conditions
    %it is usually named 'm' c.f. ''On the bvp for systems of ordinary
    %2nd order diff. equ. with a sing. of the first kind''
    
    if max(str2num(antwort_orders))==2
    %in the transformed bvp we have only separated bcs. 
    Ba=cell(nbc+2*n,6*n); %only the first nbc rows are filled
    Bb=cell(nbc+2*n,6*n); % only the last 2*n rows are filled
    else
    Ba=cell(nbc+n,2*n); 
    Bb=cell(nbc+n,2*n); 
        
    end 
    
    
    for i=1:nbc

        
        %Ba = Ba_old|Bb_old
        %        0      0  
        for j=1:n %only dirichlet werden übernommen

            Ba{i,j} = strcat(Ba_col{i}(j));

        end

        for j=1:n

           
            Ba{i,n+j} = strcat(Bb_col{i}{j})

           
        end

        if max(str2num(antwort_orders))==2
            %deal with derivative
            %setze koeff auf 0 für n+1,..2n Komponente
            for j=1:n

                Ba{i,j}=strcat(Ba{i,j},'+',Ba_col{i}(n+j),'*(-3/2)');
                Ba{i,j+2*n}=strcat('+',Ba_col{i}(n+j),'*(2)');
                Ba{i,j+4*n}=strcat('+',Ba_col{i}(n+j),'*(-1/2)');

                Ba{i,j+3*n} = '0';
                Ba{i,j+5*n} = '0';

            end

        end
        

    end %i

  

    %matrix Bb

    for i=1:n


        for j=1:n %smoothness requirementes at 1; resulting from transformation

            if j==i
                %Bb{2*n+i,j+4*n} = '1';
                 
                if max(str2num(antwort_orders))==2
                 Bb{nbc+i,j+4*n} = '1';
                else
                   Bb{nbc+i,j} = '1';  
                end 
                 
            end


        end

        for j=1:n

            if j == i
                %Bb{2*n+i,j+5*n} = '-1';
               if max(str2num(antwort_orders))==2
                Bb{nbc+i,j+5*n} = '-1';
               else 
                Bb{nbc+i,j+n} = '-1';
               end     
            end

        end
    end

    
    if max(str2num(antwort_orders))==2
    
        for i=n+1:2*n
            %setze koeff auf 0 für n+1,..2n Komponente
            %for j=1:n

            Bb{nbc+i,i-n}=strcat('+','1/2');
            Bb{nbc+i,i+ n}=strcat('-','2');
            Bb{nbc+i,i+3*n}=strcat('+','3/2');

            Bb{nbc+i,i} = strcat(Bb{i,j},'+','1/2');
            Bb{nbc+i,i+2*n} = strcat('-','2');
            Bb{nbc+i,i+4*n} = strcat('+','3/2');

            %end
        end

    end

else %antwort_endpoint ~=0
   
        
    if max(str2num(antwort_orders))==1
        Ba=cell(nbc,n); %only the first nbc rows are filled
        Bb=cell(nbc,n); % only the last 2*n rows are filled
        
        
        for i=1:nbcs
        
            for j=1:n %only dirichlet werden übernommen

                Ba{i,j} = strcat(Bb_col{i}(j));
                Bb{i,j} = strcat(Ba_col{i}(j));

            end

        
        end 

    else

        Ba=cell(nbc,3*n); %only the first nbc rows are filled
        Bb=cell(nbc,3*n); % only the last 2*n rows are filled

        for i=1:nbc



            for j=1:n %only dirichlet werden übernommen

                Ba{i,j} = strcat(Bb_col{i}(j));
                Bb{i,j} = strcat(Ba_col{i}(j));

            end

            for j=1:n
                    Ba{i,j} = strcat(Bb_col{i}(j));
                    Ba{i,j}=strcat(Ba{i,j},'+',Bb_col{i}(n+j),'*(-3/2)');
                    Ba{i,j+n}=strcat('+',Bb_col{i}(n+j),'*(2)');
                    Ba{i,j+2*n}=strcat('+',Bb_col{i}(n+j),'*(-1/2)');
            end


            

            
                for j=1:n


                    Bb{i,2*n+j} = strcat(Ba_col{i}(j));

                    Bb{i,j}=strcat('+',Ba_col{i}(n+j),'*(1/2)');
                    Bb{i,j+n}=strcat('+',Ba_col{i}(j+n),'*(-2)');
                    Bb{i,j+2*n}=strcat(Bb{i,2*n+j},'+',Ba_col{i}(j+n),'*(3/2)');


                end



        end %i


    end %first order
    
    end %antwort endpoint    

    antwort_Ba='[';

    for i=1:size(Ba,1)


        for j=1:size(Ba,2)-1


            if isempty(Ba{i,j})
                antwort_Ba=strcat(antwort_Ba,'0',',');
            else
                antwort_Ba=strcat(antwort_Ba,Ba{i,j},',');
            end

        end

        if isempty(Ba{i,size(Ba,2)})
            antwort_Ba=strcat(antwort_Ba,'0',';');
        else
            antwort_Ba=strcat(antwort_Ba,Ba{i,size(Ba,2)},';');
        end

    end

    antwort_Ba=strcat(antwort_Ba,']');
    antwort_Ba=antwort_Ba{1}

    antwort_Bb='[';

    for i=1:size(Bb,1)


        for j=1:size(Bb,2)-1

            if isempty(Bb{i,j})

                antwort_Bb=strcat(antwort_Bb,'0',',');
            else
                antwort_Bb=strcat(antwort_Bb,Bb{i,j},',');
            end
        end

        if isempty(Bb{i,size(Bb,2)})

            antwort_Bb=strcat(antwort_Bb,'0',';');
        else

            antwort_Bb=strcat(antwort_Bb,Bb{i,size(Bb,2)},';');
        end
    end

    antwort_Bb=strcat(antwort_Bb,']');
    
    if iscell(antwort_Bb)
    antwort_Bb=antwort_Bb{1};
    end 


end %inf






%if user inserts matrices, then the boundary conditions for z'(a) and z'(b)
%have to changed according to the approximation

%Call function to write input into the evpfile

if strcmp(antwort_infinite,'true')


    evpfileschreiben(filename,pfad,antwort_G,antwort_A,antwort_Ba,antwort_Bb,antwort_orders_old,N,a,b,antwort_A1,antwort_infinite,antwort_A_old,antwort_A1_old,antwort_G_old,antwort_Ba_old,antwort_Bb_old,antwort_endpoint,up_ind,low_ind);

else

    evpfileschreiben(filename,pfad,antwort_G,antwort_A,antwort_Ba,antwort_Bb,antwort_orders,N,a,b,antwort_A1,antwort_infinite,antwort_A_old,antwort_A1_old,antwort_G_old,antwort_Ba_old,antwort_Bb_old,antwort_endpoint,up_ind,low_ind);

end

evpfile = strrep(filename,'.m','');



%the function evpfileschreiben writes the input, which was read from the GUI
%in matrix form into the evpfile with the name 'evp_filename'
function ret  = evpfileschreiben(filename,pfad,antwort_G,antwort_A,antwort_Ba,antwort_Bb,antwort_orders,N,a,b,antwort_A1,antwort_infinite,antwort_A_old,antwort_A1_old,antwort_G_old,antwort_Ba_old,antwort_Bb_old,antwort_endpoint,up_ind,low_ind)

if strcmp(antwort_infinite,'true')
    n=length(str2num(antwort_orders));
   
    if antwort_endpoint==0
    n=2*n;
    
    end
else
    n=length(str2num(antwort_orders));
end
n=num2str(n);
if length(pfad) >0
    datei=fopen(strcat(pfad,filename),'w');
else
    datei=fopen(filename,'w');
end
help = strrep(filename,'.m','');



fprintf(datei,'function [ret] = %s(nr,t)\n',help);
fprintf(datei,'%%Kontrollnummer43753976430976655145\n');

fprintf(datei,'\n switch nr \n');
mainGUIhandle=bvpsuite;
mainGUIdata  = guidata(mainGUIhandle);



if get(mainGUIdata.evp_mm,'Value')==1
    antwort_EVP_mm='true';
else
    antwort_EVP_mm='false';
end

fprintf(datei,'\n    case ''EVP_mm''\n          ret=%s;\n',char(antwort_EVP_mm));
fprintf(datei,'\n    case ''G''\n          ret=%s;\n',antwort_G);
fprintf(datei,'\n    case ''A0''\n          ret=%s;\n',antwort_A);
fprintf(datei,'\n    case ''A1''\n          ret=%s;\n',antwort_A1);
fprintf(datei,'\n    case ''Ba''\n          ret=%s;\n',antwort_Ba);
fprintf(datei,'\n    case ''Bb''\n          ret=%s;\n',antwort_Bb);
fprintf(datei,'\n    case ''a''\n          ret=%i;\n',a);
    fprintf(datei,'\n    case ''b''\n          ret=%i;\n',b);
    fprintf(datei,'\n    case ''N''\n          ret=%i;\n',N);
    if strcmp(antwort_infinite,'true')
        fprintf(datei,'\n    case ''antwort_endpoint''\n          ret=%i;\n',antwort_endpoint);
    else
        fprintf(datei,'\n    case ''antwort_endpoint''\n          ret=[];\n');
    end
% else
%     fprintf(datei,'\n case ''a''\n  ret=%s;\n',char(a));
%     fprintf(datei,'\n case ''b''\n  ret=%s;\n',char(b));
%     fprintf(datei,'\n case ''N''\n  ret=%s;\n',char(N));
% 
% end

fprintf(datei,'\n    case ''infinite''\n          ret=%s;\n',char(antwort_infinite));
fprintf(datei,'\n    case ''n''\n          ret=%s;\n',n);
fprintf(datei,'\n    case ''orders''\n          ret=%s;\n',antwort_orders);
fprintf(datei,'end \n')
fprintf(datei,'%%Values read by bvpsuite GUI:\n');
fprintf(datei,'%%EVP_mm\n%%%s%%#EVP_mm\n',char(antwort_EVP_mm));
fprintf(datei,'%%Values read by evpsolve GUI:\n');
fprintf(datei,'%%infinite\n%%%s%%#infinite\n',char(antwort_infinite));
fprintf(datei,'%%orders\n%%%s%%#orders\n',antwort_orders);
%fprintf(datei,'%%discret\n%%%s%%#discret\n',char(antwort_discret));

if strcmp(antwort_infinite,'true')
    if antwort_endpoint ~=0
        fprintf(datei,'%%a\n%%%s%%#a\n',num2str(antwort_endpoint));
    else
        fprintf(datei,'%%a\n%%%s%%#a\n',num2str(a));
    end
    
end
fprintf(datei,'%%b\n%%%s%%#b\n',num2str(b));
fprintf(datei,'%%N\n%%%s%%#N\n',num2str(N));
fprintf(datei,'%%A0\n%%%s%%#A0\n',antwort_A_old);
fprintf(datei,'%%A1\n%%%s%%#A1\n',antwort_A1_old);
fprintf(datei,'%%G\n%%%s%%#G\n',antwort_G_old);
fprintf(datei,'%%Ba\n%%%s%%#Ba\n',antwort_Ba_old);
fprintf(datei,'%%Bb\n%%%s%%#Bb\n',antwort_Bb_old);
fprintf(datei,'%%up_ind\n%%%s%%#up_ind\n',num2str(up_ind));
fprintf(datei,'%%low_ind\n%%%s%%#low_ind\n',num2str(low_ind));
fprintf(1,'The file %s was written!\n',help);
fprintf(1,'Check the inputs in the file!\n');


evpsolve
helpdlg('The file has been written!','SUCCESS');


%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton5.
function pushbuttonSolve_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


low_ind=get(handles.low_ind,'String');
up_ind=get(handles.up_ind,'String');

low_ind=str2num(low_ind);
up_ind=str2num(up_ind);




mainGUIhandle=bvpsuite;
mainGUIdata  = guidata(mainGUIhandle);


filename=get(mainGUIdata.edit1,'String');
pfad=get(mainGUIdata.pfad,'String')



%filename=strcat('evp_',filename)
evpfile = strrep(filename,'.m','')



antwort_infinite=get(handles.radiobuttonInfinite,'Value');
antwort_orders= get(handles.editOrders,'String');


aktuellesverzeichnis=cd;


 if  max(str2num(antwort_orders)) == 1
 [evhelp,lambdas,A,B] = boxscheme(pfad,evpfile,5);
 else
[evhelp,lambdas,A,B] = approx2ndorder(pfad,evpfile,5);
 end

 
 cd(aktuellesverzeichnis);
 
t =sym('t');
%if det(A-t*B) == 0


%    errordlg('The resulting algebraic eigenvalue problem is singular!')

% else
orders=str2num(antwort_orders);

%Calls for the functions boxscheme/approx2ndorder
%to set up and solve the algebraic eigenvalue problem
%see the corresponding m-files for detailed information


if orders(1) == 1
    [evhelp,lambdas,A,B,xhelp] = boxscheme(pfad,evpfile);

else
    [evhelp,lambdas,A,B,xhelp] = approx2ndorder(pfad,evpfile);
end

 cd(aktuellesverzeichnis);

 
if get(handles.printev,'Value')==1 


   if isempty(low_ind) | isempty(up_ind) | (isempty(low_ind) & isempty(up_ind))    
      err_evp('evp1');
      return;
   end    
%prints all the eigenvalues in a sorted manner
fprintf('\n\\  \\begin{tabular}{|c|}\n');
fprintf('\\hline\n');
fprintf(' $\\lambda $\\\\ \n');
fprintf('\\hline\n');
for i=low_ind:up_ind
    fprintf ( 1, '%s\n', num2str(lambdas(i)) );
end
fprintf('\\hline\n');
fprintf('\\end{tabular} \n\n');


end


figure(evpsolve)


if get(handles.plotef,'Value')==1 
 
   if isempty(low_ind) | isempty(up_ind) | (isempty(low_ind) & isempty(up_ind))    
      err_evp('evp1');
      return;
   end    
    
  figure
 
  
   diffindex=up_ind-low_ind+1;
 
   for i=1:diffindex
       
         for j=1:length(orders)    
           
                         
            subplot(diffindex*length(orders),1,(i-1)*length(orders)+j);  
            plot(xhelp,evhelp(:,(low_ind-1)*length(orders) + (i-1)*length(orders) +j));
            
            
            
         end 
    
    title(subplot(diffindex*length(orders),1,1+(i-1)*length(orders)),['Eigenfunction for eigenvalue \lambda = ',num2str(lambdas(low_ind+i-1))]);
  end     
    

  
end     


helpdlg('The calculation was successful. The eigenvalues and eigenfunctions are saved in evplast.mat. The eigenfunctions correspond to the column vectors of the matrix ''ev''. The eigenvalues are listed in ''lambdas''. ''x'' denotes the discrete vector of the independent variable. You can use them now in the workspace (For further information see "Workspace" in the Matlab help!).','SUCCESS');



% --- Executes on button press in pushbutton6.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mainGUIhandle=bvpsuite;
mainGUIdata  = guidata(mainGUIhandle);

set(mainGUIdata.speichern,'Enable','on');


close(evpsolve)



% --- Executes on button press in pushbutton7.
function pushbuttonDelete_Callback(hObject, eventdata, handles)

% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.editOrders,'String','');
set(handles.editA1,'String','');
set(handles.editA0,'String','');
set(handles.editG,'String','');
set(handles.editB0,'String','');
set(handles.editB1,'String','');
%set(handles.editEndpoint,'String','');
set(handles.edita,'String','');
set(handles.editb,'String','');
set(handles.editN,'String','');
set(handles.radiobuttonInfinite,'Value',0);
set(handles.low_ind,'String','');
set(handles.up_ind,'String','');
set(handles.plotef,'Value',0);



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbuttonMesh_Callback(hObject, eventdata, handles)
helpdlg('Insert the mesh in the form [a,b,N] whereas N denotes the number of intervals on [a,b].','Info');
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





% --- Executes on button press in pushbuttonLoad.
function pushbuttonLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


mainGUIhandle=bvpsuite;
mainGUIdata  = guidata(mainGUIhandle);


filename=get(mainGUIdata.edit1,'String');
pfad=get(mainGUIdata.pfad,'String');
%filename = strcat('evp_',filename);
%evpfile = strrep(filename,'.m','')
datei=fopen(strcat(pfad,filename),'r');
if (datei~=-1)

    %C:wenn man ein bestehendes file neu laedt

    %C: speichere inhalt von datei in inhaltdatei
    inhaltdatei=fread(datei);
    inhaltdatei=char(inhaltdatei');

    help=inhaltdatei(strfind(inhaltdatei,'%infinite'):strfind(inhaltdatei,'%#infinite')-1);
    help=strrep(help,'%infinite','');
    help(3:length(help));
    help=help(3:length(help));
    infinite=help;

    help=inhaltdatei(strfind(inhaltdatei,'%orders'):strfind(inhaltdatei,'%#orders')-1);
    help=strrep(help,'%orders','');
    help(3:length(help));
    help=help(3:length(help));
    orders=help;

    help=inhaltdatei(strfind(inhaltdatei,'%A0'):strfind(inhaltdatei,'%#A0')-1);
    help=strrep(help,'%A0','');
    help(3:length(help));
    help=help(3:length(help));
    A0=help;

    help=inhaltdatei(strfind(inhaltdatei,'%A1'):strfind(inhaltdatei,'%#A1')-1);
    help=strrep(help,'%A1','');
    help(3:length(help));
    help=help(3:length(help));
    A1=help;

    help=inhaltdatei(strfind(inhaltdatei,'%G'):strfind(inhaltdatei,'%#G')-1);
    help=strrep(help,'%G','');
    help(3:length(help));
    help=help(3:length(help));
    G=help;

    help=inhaltdatei(strfind(inhaltdatei,'%Ba'):strfind(inhaltdatei,'%#Ba')-1);
    help=strrep(help,'%Ba','');
    help(3:length(help));
    help=help(3:length(help));
    Ba=help;

    help=inhaltdatei(strfind(inhaltdatei,'%Bb'):strfind(inhaltdatei,'%#Bb')-1);
    help=strrep(help,'%Bb','');
    help(3:length(help));
    help=help(3:length(help));
    Bb=help;

    help=inhaltdatei(strfind(inhaltdatei,'%a'):strfind(inhaltdatei,'%#a')-1);
    help=strrep(help,'%a','');
    help=help(3:length(help));
    a=help;
    
       
    help=inhaltdatei(strfind(inhaltdatei,'%b'):strfind(inhaltdatei,'%#b')-1);
    help=strrep(help,'%b','');
    help=help(3:length(help));
    b=help;
    
    help=inhaltdatei(strfind(inhaltdatei,'%N'):strfind(inhaltdatei,'%#N')-1);
    help=strrep(help,'%N','');
    help=help(3:length(help));
    N=help;
    
    help=inhaltdatei(strfind(inhaltdatei,'%up_ind'):strfind(inhaltdatei,'%#up_ind')-1);
    help=strrep(help,'%up_ind','');
    help=help(3:length(help));
    up_ind=help;
    
    help=inhaltdatei(strfind(inhaltdatei,'%low_ind'):strfind(inhaltdatei,'%#low_ind')-1);
    help=strrep(help,'%low_ind','');
    help=help(3:length(help));
    low_ind=help;
    
else    
    errordlg('The corresponding evpfile does not exist yet! Check if your problem is an eigenvalue problem! ','Error');
  
    return;



end

set(handles.editOrders,'String',orders);

set(handles.edita,'String',a)


set(handles.editN,'String',N)



    set(handles.editA0,'String',A0);
    set(handles.editA1,'String',A1);
    set(handles.editG,'String',G);
    set(handles.editB0,'String',Ba);
    set(handles.editB1,'String',Bb);
    set(handles.up_ind,'String',up_ind);
    set(handles.low_ind,'String',low_ind);



checkInfinite= length(strfind(infinite,'true'));

if checkInfinite ~= 0
    set(handles.radiobuttonInfinite,'Value',1);
     set(handles.editb,'String','infinity');
     set(handles.editb,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   
else
    set(handles.radiobuttonInfinite,'Value',0);
    set(handles.editb,'String',b);
end


 set(mainGUIdata.speichern,'Enable','off');
  

figure(evpsolve)





function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editA1_Callback(hObject, eventdata, handles)
% hObject    handle to editA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% Hints: get(hObject,'String') returns contents of editA1 as text
%        str2double(get(hObject,'String')) returns contents of editA1 as a double


% --- Executes during object creation, after setting all properties.
function editA1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Your system has either to be of second order, -y'''' + A1y'' + A0y = lambdaGy, or of first order, -y'' +A0y = lambdaGy. Insert the matrices A1, A0 and G in the corresponding fields. If A1 or A0 is empty, insert a zero matrix of the corresponding size. See the manual within the diploma thesis of C. Simon for more information.','help');



% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('The boundary conditions have to be linear, i.e. Bay(a) +Bby(b)=0. Insert the constant matrices Ba and Bb. If the differential operator is of order 2, then Ba and Bb are of dimension (2n)x(2n), where n is the number of equations.','help')




% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2





function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a
%        double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
helpdlg('If your problem is posed on a semi-infinite interval tick the check box to the left.','help');

% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)








% --- Executes on button press in pushbuttonLoad.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Insert the number of subintervals used for the discretization. In case of a semi-infinite interval ''N'' defines the number of subintervals on [0,1]. Therefore the solution on the truncated interval is provided on 2N-1 points.','help');




function editb_Callback(hObject, eventdata, handles)
% hObject    handle to editb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editb as text
%        str2double(get(hObject,'String')) returns contents of editb as a double


% --- Executes during object creation, after setting all properties.
function editb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editN_Callback(hObject, eventdata, handles)
% hObject    handle to editN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editN as text
%        str2double(get(hObject,'String')) returns contents of editN as a double


% --- Executes during object creation, after setting all properties.
function editN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edita_Callback(hObject, eventdata, handles)
% hObject    handle to edita (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edita as text
%        str2double(get(hObject,'String')) returns contents of edita as a double




% --- Executes during object creation, after setting all properties.
function edita_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edita (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function low_ind_Callback(hObject, eventdata, handles)
% hObject    handle to low_ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_ind as text
%        str2double(get(hObject,'String')) returns contents of low_ind as a double


% --- Executes during object creation, after setting all properties.
function low_ind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function up_ind_Callback(hObject, eventdata, handles)
% hObject    handle to up_ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of up_ind as text
%        str2double(get(hObject,'String')) returns contents of up_ind as a double


% --- Executes during object creation, after setting all properties.
function up_ind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to up_ind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('The eigenvalues with indices between the two inserted integers are printed in the Matlab command window. They are sorted in ascending order. If the spectrum contains complex eigenvalues, whose imaginary part is unequal 0, they are listed at the end. Complex numbers are sorted by their absolute values and matches are further sorted by their angles..','help')



% --- Executes on button press in plotef.
function plotef_Callback(hObject, eventdata, handles)
% hObject    handle to plotef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotef
%if get(handles.plotef,'Value')==1
    
    
%     antwort_orders= get(handles.editOrders,'String');
%     orders=str2num(char(antwort_orders));
%     datei=fopen('lastevp.mat');
%     if datei==-1
%         error('The file "evplast.mat" does not exist, run a calculation to create it ');
%         return;
%     end
%     
%      fclose(datei);
%      load lastevp;
%      
%      low_ind=get(handles.low_ind,'String');
%      low_ind=str2num(low_ind);
%      up_ind=get(handles.up_ind,'String');
%      up_ind=str2num(up_ind);
%     
%      if isempty(low_ind) | isempty(up_ind) | (isempty(low_ind) & isempty(up_ind))    
%       err_evp('evp1');
%       return;
%      end    
%      
%      
%     figure
%     diffindex=up_ind-low_ind+1;
%  
%   
%     
%    for i=1:diffindex
%        
%          for j=1:length(orders)    
%            
%                          
%             subplot(diffindex*length(orders),1,(i-1)*length(orders)+j);  
%             plot(x,ev(:,(low_ind-1)*length(orders) + (i-1)*length(orders) +j));
%             
%             
%             
%          end 
%     
%     title(subplot(diffindex*length(orders),1,1+(i-1)*length(orders)),['Eigenfunction for eigenvalue \lambda = ',num2str(lambdas(low_ind+i-1))]);
%   end     
    
  
       
       
       
   
  
   

%end 

% --- Executes on button press in printev.
function printev_Callback(hObject, eventdata, handles)
% hObject    handle to printev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of printev


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
helpdlg('Insert the left endpoint ''a'' and the right endpoint ''b'' of the interval on which the problem is posed. If your problem is posed on a semi-infinite interval, leave ''b'' empty.','help')
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Tick this checkbox if you want to have a plot of the corresponding eigenfunctions.','help');



% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datei=fopen('lastevp.mat');
    if datei==-1
       errordlg('The file "lastevp.mat" does not exist. Run a calculation to create it!');
        return;
    end
    fclose(datei);
load lastevp.mat;
assignin('base','x',x);
assignin('base','ev',ev);
assignin('base','lambdas',lambdas);
open x;
open ev;
open lambdas;


