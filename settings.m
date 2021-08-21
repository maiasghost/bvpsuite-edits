function varargout = einstellungen(varargin)
% EINSTELLUNGEN M-file for einstellungen.fig
%      EINSTELLUNGEN, by itself, creates a new EINSTELLUNGEN or raises the existing
%      singleton*.
%
%      H = EINSTELLUNGEN returns the handle to a new EINSTELLUNGEN or the handle to
%      the existing singleton*.
%
%      EINSTELLUNGEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EINSTELLUNGEN.M with the given input arguments.
%
%      EINSTELLUNGEN('Property','Value',...) creates a new EINSTELLUNGEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before einstellungen_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to einstellungen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help einstellungen

% Last Modified by GUIDE v2.5 18-Aug-2009 22:18:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @einstellungen_OpeningFcn, ...
                   'gui_OutputFcn',  @einstellungen_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before einstellungen is made visible.
function einstellungen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to einstellungen (see VARARGIN)

% Choose default command line output for einstellungen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes einstellungen wait for user response (see UIRESUME)
% uiwait(handles.figure1);
load options
set(handles.abstol,'String',num2str(abstol));
set(handles.reltol,'String',num2str(reltol));
set(handles.maxiter,'String',num2str(maxiter));
set(handles.maxfunevals,'String',num2str(maxfunevals));
set(handles.updatejacfactor,'String',num2str(updatejacfactor));
set(handles.lambdamin,'String',num2str(lambdamin));
set(handles.switchtoffnfactor,'String',num2str(switchtoffnfactor));
set(handles.ausgabe,'String',num2str(ausgabe));
set(handles.TRM,'Value',TRM);
set(handles.abstolgitter,'String',num2str(abstolgitter));
set(handles.reltolgitter,'String',num2str(reltolgitter));
set(handles.K,'String',num2str(K));
set(handles.wiederholung,'String',num2str(wiederholung));
set(handles.feinesgitter,'Value',feinesgitter);
set(handles.n0,'String',num2str(n0));


% --- Outputs from this function are returned to the command line.
function varargout = einstellungen_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in speichern.
function speichern_Callback(hObject, eventdata, handles)
abstol=str2num(get(handles.abstol,'String'));
reltol=str2num(get(handles.reltol,'String'));
maxiter=str2num(get(handles.maxiter,'String'));
maxfunevals=str2num(get(handles.maxfunevals,'String'));
updatejacfactor=str2num(get(handles.updatejacfactor,'String'));
lambdamin=str2num(get(handles.lambdamin,'String'));
switchtoffnfactor=str2num(get(handles.switchtoffnfactor,'String'));
TRM=get(handles.TRM,'Value');
abstolgitter=str2num(get(handles.abstolgitter,'String'));
reltolgitter=str2num(get(handles.reltolgitter,'String'));
K=str2num(get(handles.K,'String'));
n0=str2num(get(handles.n0,'String'));
wiederholung=str2num(get(handles.wiederholung,'String'));
default_zf_opt = optimset('Display','off','MaxIter',maxiter,'MaxFunEvals',maxfunevals);
defaultopt = sbvpset('AbsTol',abstol,'Basis','RungeKutta','CheckJac',0,'ColPts','Equidistant','Degree',4,'DegreeSelect','auto','Display',1,'fVectorized',0,'IntMaxMinRatio',10,'JacVectorized',0,'MaxMeshPts',10000,'OutputFcn','','OutputSel',[],'OutputTrace',1,'RelTol',reltol,'ZfOpt',default_zf_opt);
bvpopt = sbvpset;
bvpopt.ZfOpt = optimset(defaultopt.ZfOpt,bvpopt.ZfOpt);
bvpopt = sbvpset(defaultopt, bvpopt);
bvpopt.Log = 0;
bvpopt.Private.UpdateJacFactor = updatejacfactor;
bvpopt.Private.LambdaMin = lambdamin;
bvpopt.Private.SwitchToFFNFactor = switchtoffnfactor;
bvpopt.TRM = TRM;
ausgabe=str2num(get(handles.ausgabe,'String'));
feinesgitter=get(handles.feinesgitter,'Value');

if length(strfind(char(version),'R14'))==0
    save options bvpopt abstol reltol maxiter maxfunevals updatejacfactor lambdamin switchtoffnfactor ausgabe TRM abstolgitter reltolgitter K wiederholung feinesgitter n0;
else
    save options -v6 bvpopt abstol reltol maxiter maxfunevals updatejacfactor lambdamin switchtoffnfactor ausgabe TRM abstolgitter reltolgitter K wiederholung feinesgitter n0; 
end

load options.mat
%open abstol

% get the main_gui handle (access to the gui)
mainGUIhandle = bvpsuite ; 
% get the data from the gui (all handles inside gui_main)
mainGUIdata  = guidata(mainGUIhandle);


 


if get(mainGUIdata.automatic,'Value')==0
   set(mainGUIdata.edit13,'BackgroundColor','white');
   set(mainGUIdata.edit13,'Enable','on');
   %set(handles.n0,'String','');

elseif get(mainGUIdata.automatic,'Value')==1
   set(mainGUIdata.edit13,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   set(mainGUIdata.edit13,'Enable','off');    

  load options.mat 
  
  %open abstol
  
   
  if  abstolgitter <= 10^-6                             

    set(mainGUIdata.edit13, 'String', 8);
  
  elseif abstolgitter > 10^-6 && abstolgitter <= 10^-4
      
    set(mainGUIdata.edit13, 'String', 6);
  
  elseif abstolgitter > 10^-4 && abstolgitter <= 10^-2
      
     set(mainGUIdata.edit13, 'String', 4);
  else abstolgitter > 10^-2 
  
     set(mainGUIdata.edit13, 'String', 2);
  end     
    
  
 
    
end

% save changed data back into main_gui
%this line updates the data of the Main Gui
guidata(bvpsuite, mainGUIdata);
close(settings)



if get(mainGUIdata.automatic,'Value')==1
    helpdlg('Press ''Save'' to save your changes.','help'); 

end
% --- Executes during object creation, after setting all properties.
function abstol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to abstol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function abstol_Callback(hObject, eventdata, handles)
% hObject    handle to abstol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of abstol as text
%        str2double(get(hObject,'String')) returns contents of abstol as a double


% --- Executes during object creation, after setting all properties.
function reltol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reltol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function reltol_Callback(hObject, eventdata, handles)
% hObject    handle to reltol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reltol as text
%        str2double(get(hObject,'String')) returns contents of reltol as a double


% --- Executes on button press in standard.
function standard_Callback(hObject, eventdata, handles)
set(handles.abstol,'String','1e-10');
set(handles.reltol,'String','1e-10');
set(handles.maxiter,'String','90000000');
set(handles.maxfunevals,'String','90000000');
set(handles.updatejacfactor,'String','0.5');
set(handles.lambdamin,'String','0.001');
set(handles.switchtoffnfactor,'String','0.5');
set(handles.ausgabe,'String','1');
set(handles.TRM,'Value',0);
set(handles.feinesgitter,'Value',0);
set(handles.abstolgitter,'String',1e-6);
set(handles.reltolgitter,'String',1e-3);
set(handles.K,'String',10000);
set(handles.n0,'String',50);
set(handles.wiederholung,'String',15);

% --- Executes on button press in abbrechen.
function abbrechen_Callback(hObject, eventdata, handles)
close;


% --- Executes during object creation, after setting all properties.
function maxiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxiter_Callback(hObject, eventdata, handles)
% hObject    handle to maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxiter as text
%        str2double(get(hObject,'String')) returns contents of maxiter as a double


% --- Executes during object creation, after setting all properties.
function maxfunevals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxfunevals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxfunevals_Callback(hObject, eventdata, handles)
% hObject    handle to maxfunevals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxfunevals as text
%        str2double(get(hObject,'String')) returns contents of maxfunevals as a double


% --- Executes during object creation, after setting all properties.
function updatejacfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to updatejacfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function updatejacfactor_Callback(hObject, eventdata, handles)
% hObject    handle to updatejacfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of updatejacfactor as text
%        str2double(get(hObject,'String')) returns contents of updatejacfactor as a double


% --- Executes during object creation, after setting all properties.
function lambdamin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambdamin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function lambdamin_Callback(hObject, eventdata, handles)
% hObject    handle to lambdamin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambdamin as text
%        str2double(get(hObject,'String')) returns contents of lambdamin as a double


% --- Executes during object creation, after setting all properties.
function switchtoffnfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to switchtoffnfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function switchtoffnfactor_Callback(hObject, eventdata, handles)
% hObject    handle to switchtoffnfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of switchtoffnfactor as text
%        str2double(get(hObject,'String')) returns contents of switchtoffnfactor as a double


% --- Executes during object creation, after setting all properties.
function ausgabe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ausgabe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ausgabe_Callback(hObject, eventdata, handles)
% hObject    handle to ausgabe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ausgabe as text
%        str2double(get(hObject,'String')) returns contents of ausgabe as a double


% --- Executes on button press in TRM.
function TRM_Callback(hObject, eventdata, handles)
% hObject    handle to TRM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TRM




% --- Executes on button press in fehler.
function fehler_Callback(hObject, eventdata, handles)
% hObject    handle to fehler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fehler





function abstolgitter_Callback(hObject, eventdata, handles)
% hObject    handle to abstolgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of abstolgitter as text
%        str2double(get(hObject,'String')) returns contents of abstolgitter as a double


% --- Executes during object creation, after setting all properties.
function abstolgitter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to abstolgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function K_Callback(hObject, eventdata, handles)
% hObject    handle to K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K as text
%        str2double(get(hObject,'String')) returns contents of K as a double


% --- Executes during object creation, after setting all properties.
function K_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function reltolgitter_Callback(hObject, eventdata, handles)
% hObject    handle to reltolgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reltolgitter as text
%        str2double(get(hObject,'String')) returns contents of reltolgitter as a double


% --- Executes during object creation, after setting all properties.
function reltolgitter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reltolgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in gittersteuerung.
function gittersteuerung_Callback(hObject, eventdata, handles)
% hObject    handle to gittersteuerung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gittersteuerung





function wiederholung_Callback(hObject, eventdata, handles)
% hObject    handle to wiederholung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wiederholung as text
%        str2double(get(hObject,'String')) returns contents of wiederholung as a double


% --- Executes during object creation, after setting all properties.
function wiederholung_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wiederholung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in feinesgitter.
function feinesgitter_Callback(hObject, eventdata, handles)
% hObject    handle to feinesgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of feinesgitter





function N0_Callback(hObject, eventdata, handles)
% hObject    handle to N0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N0 as text
%        str2double(get(hObject,'String')) returns contents of N0 as a double


% --- Executes during object creation, after setting all properties.
function N0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function n0_Callback(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n0 as text
%        str2double(get(hObject,'String')) returns contents of n0 as a double


% --- Executes during object creation, after setting all properties.
function n0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


