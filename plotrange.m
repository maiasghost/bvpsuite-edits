function varargout = plotrange(varargin)
% PLOTRANGE M-file for plotrange.fig
%      PLOTRANGE, by itself, creates a new PLOTRANGE or raises the existing
%      singleton*.
%
%      H = PLOTRANGE returns the handle to a new PLOTRANGE or the handle to
%      the existing singleton*.
%
%      PLOTRANGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTRANGE.M with the given input arguments.
%
%      PLOTRANGE('Property','Value',...) creates a new PLOTRANGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plotrange_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plotrange_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plotrange

% Last Modified by GUIDE v2.5 08-May-2009 10:31:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotrange_OpeningFcn, ...
                   'gui_OutputFcn',  @plotrange_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before plotrange is made visible.
function plotrange_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plotrange (see VARARGIN)

% Choose default command line output for plotrange
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plotrange wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plotrange_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function leftendpoint_Callback(hObject, eventdata, handles)
% hObject    handle to leftendpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of leftendpoint as text
%        str2double(get(hObject,'String')) returns contents of leftendpoint as a double


% --- Executes during object creation, after setting all properties.
function leftendpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to leftendpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rightendpoint_Callback(hObject, eventdata, handles)
% hObject    handle to rightendpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rightendpoint as text
%        str2double(get(hObject,'String')) returns contents of rightendpoint as a double


% --- Executes during object creation, after setting all properties.
function rightendpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rightendpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load last.mat

%plotpoly(polynomials,x1,length(valx1(:,1)));
n=length(valx1tau(:,1));

main_gui_handle=bvpsuite;
mainGUIdata  = guidata(main_gui_handle);



pfad=get(mainGUIdata.pfad,'String');

aktuellesverzeichnis=cd;
arbeitsverzeichnis=pfad;
cd(arbeitsverzeichnis);


filename=get(mainGUIdata.edit1,'String');
bvpfile=strrep(filename,'.m','');



%bvpfile=strcat(pfad,bvpfile)

left=str2num(get(handles.leftendpoint,'String'));
right=str2num(get(handles.rightendpoint,'String'));

if feval(bvpfile,'Infinite')==1 &  feval(bvpfile,'EVP')==1 



    ep = feval(bvpfile,'Endpoint');


    if ep == 0

        cd(aktuellesverzeichnis);
        for i=1:n/2

            
            
            [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);



        end


        plot_results(tau_infinite,sol_infinite,n/2-2,1,1,left,right)
        

    else

         cd(aktuellesverzeichnis);
         
     
         
         
        for i=1:n

            [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

        end


         plot_results(tau_infinite,sol_infinite,n-2,1,1,left,right)
        
        

    end





elseif feval(bvpfile,'Infinite')==1 &  feval(bvpfile,'EVP')==0

    ep = feval(bvpfile,'Endpoint');

    if ep==0

         cd(aktuellesverzeichnis);
         
        for i=1:n/2


            [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);

        end
        
         plot_results(tau_infinite,sol_infinite,n/2,0,1,left,right)


    else

         cd(aktuellesverzeichnis); 
         
        for i=1:n

            [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

        end
        
        
         plot_results(tau_infinite,sol_infinite,n,0,1,left,right)



    end

end

close plotrange

