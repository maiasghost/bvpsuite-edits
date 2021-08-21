function varargout = pathfollowinggui(varargin)
% PATHFOLLOWINGGUI M-file for pathfollowinggui.fig
%      PATHFOLLOWINGGUI, by itself, creates a new PATHFOLLOWINGGUI or 
%      raises the existing
%      singleton*.
%
%      H = PATHFOLLOWINGGUI returns the handle to a new PATHFOLLOWINGGUI or the handle to
%      the existing singleton*.
%
%      PATHFOLLOWINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PATHFOLLOWINGGUI.M with the given input arguments.
%
%      PATHFOLLOWINGGUI('Property','Value',...) creates a new PATHFOLLOWINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pathfollowinggui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pathfollowinggui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help pathfollowinggui

% Last Modified by GUIDE v2.5 09-Jul-2010 20:37:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pathfollowinggui_OpeningFcn, ...
                   'gui_OutputFcn',  @pathfollowinggui_OutputFcn, ...
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


% --- Executes just before pathfollowinggui is made visible.
function pathfollowinggui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pathfollowinggui (see VARARGIN)


% Choose default command line output for pathfollowinggui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pathfollowinggui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if length(strfind(char(cd),'\'))>0
    set(handles.pfad,'String',strcat(char(cd),'\')); %f�r Windows
else
    set(handles.pfad,'String',strcat(char(cd),'/')); %f�r Linux/Unix
end
set(handles.name,'String','last.mat');

if length(strfind(char(cd),'\'))>0
    set(handles.speicherpfad,'String',strcat(char(cd),'\')); %f�r Windows
else
    set(handles.speicherpfad,'String',strcat(char(cd),'/')); %f�r Linux/Unix
end

set(handles.speichername,'String','.mat');
set(handles.schrittanzahl,'String','1');
set(handles.schrittweite,'String','0.1');
set(handles.startindex,'String','1');
set(handles.meshadaptation,'Value',1);
set(handles.pfaddaten,'String','[1 0]');
load help1;
%set(handles.zeichnen,'Value',zeichnen);



% --- Outputs from this function are returned to the command line.
function varargout = pathfollowinggui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function pfad_Callback(hObject, eventdata, handles)
% hObject    handle to pfad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pfad as text
%        str2double(get(hObject,'String')) returns contents of pfad as a
% double


% --- Executes during object creation, after setting all properties.
function pfad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pfad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name as text
%        str2double(get(hObject,'String')) returns contents of name as a
% double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in berechnen.
function berechnen_Callback(hObject, eventdata, handles)

pfad=char(get(handles.pfad,'String'));
matfile=char(get(handles.name,'String'));
speicherpfad=char(get(handles.speicherpfad,'String'));
schrittanzahl=str2num(get(handles.schrittanzahl,'String'));
speichername=char(get(handles.speichername,'String'));
speichername2=speichername;
if length(strfind(speichername,'.mat'))>0
    speichername2=speichername(1:strfind(speichername,'.mat')-1);
end
schrittweite=str2num(get(handles.schrittweite,'String'));
startindex=str2num(get(handles.startindex,'String'));
hitpoint=str2num(get(handles.hitpoint,'String'));
pfaddaten=str2num(get(handles.pfaddaten,'String'));
meshadaptation=get(handles.meshadaptation,'Value');
%Wenn der Startindex gr��er als 1 gesetzt ist, m�chte man mit einer bereits
%per Pathfollowing berechneten L�sung beginnen, diese holt sich das
%Programm von alleine
if startindex>1
    matfile=strcat(speichername2,num2str(startindex));
    load(strcat(strcat(speicherpfad,matfile),'.mat'));
else
    load(strcat(pfad,matfile));
end
load help1;
load options;
zaehler=[];
zeichnen=0;

pathfollowing(bvpfile,zeichnen,x1,coeff,schrittweite,schrittanzahl,bvpopt,ausgabe,speicherpfad,startindex,speichername,pfaddaten,meshadaptation,feinesgitter,hitpoint);





function speicherpfad_Callback(hObject, eventdata, handles)
% hObject    handle to speicherpfad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speicherpfad as text
%        str2double(get(hObject,'String')) returns contents of speicherpfad
% as a double


% --- Executes during object creation, after setting all properties.
function speicherpfad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speicherpfad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in suchen.
function suchen_Callback(hObject, eventdata, handles)

aktuellesverzeichnis=cd;
arbeitsverzeichnis=get(handles.pfad,'String');
cd(arbeitsverzeichnis);

[FileName,PathName,index] = uigetfile('*.mat','Choose your mat-file');

if index>0
    set(handles.name,'String',FileName);
    set(handles.pfad,'String',PathName);
end
cd(aktuellesverzeichnis);



% --- Executes on button press in suchen2.
function suchen2_Callback(hObject, eventdata, handles)



aktuellesverzeichnis=cd;
arbeitsverzeichnis=get(handles.pfad,'String');
cd(arbeitsverzeichnis);

verzeichnis = uigetdir(arbeitsverzeichnis,'Where do you want to save your data?');
if length(strfind(char(cd),'\'))>0
    set(handles.speicherpfad,'String',strcat(verzeichnis,'\')); %f�r Windows
else
    set(handles.speicherpfad,'String',strcat(verzeichnis,'/')); %f�r
Linux/Unix
end

cd(aktuellesverzeichnis);




function schrittanzahl_Callback(hObject, eventdata, handles)
% hObject    handle to schrittanzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of schrittanzahl as text
%        str2double(get(hObject,'String')) returns contents of schrittanzahl as a double


% --- Executes during object creation, after setting all properties.
function schrittanzahl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to schrittanzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function speichername_Callback(hObject, eventdata, handles)
% hObject    handle to speichername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of speichername as text
%        str2double(get(hObject,'String')) returns contents of speichername as a double


% --- Executes during object creation, after setting all properties.
function speichername_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speichername (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function schrittweite_Callback(hObject, eventdata, handles)
% hObject    handle to schrittweite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of schrittweite as text
%        str2double(get(hObject,'String')) returns contents of schrittweite as a double


% --- Executes during object creation, after setting all properties.
function schrittweite_CreateFcn(hObject, eventdata, handles)
% hObject    handle to schrittweite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in pfaddarstellen.
function pfaddarstellen_Callback(hObject, eventdata, handles)

speicherpfad=get(handles.speicherpfad,'String');
speichername=get(handles.speichername,'String');
load(strcat(speicherpfad,speichername));
figure;
plot(parametervalue,pathdata,'color','black');
xlabel('p1');
ylabel('functional');




function startindex_Callback(hObject, eventdata, handles)
% hObject    handle to startindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startindex as text
%        str2double(get(hObject,'String')) returns contents of startindex as a double


% --- Executes during object creation, after setting all properties.
function startindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function pfadddaten_Callback(hObject, eventdata, handles)
% hObject    handle to pfadddaten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pfadddaten as text
%        str2double(get(hObject,'String')) returns contents of pfadddaten as a double


% --- Executes during object creation, after setting all properties.
function pfadddaten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pfadddaten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in meshadaptation.
function meshadaptation_Callback(hObject, eventdata, handles)
% hObject    handle to meshadaptation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meshadaptation




% --- Executes on button press in zeichnen.
function zeichnen_Callback(hObject, eventdata, handles)
% hObject    handle to zeichnen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zeichnen

%zeichnen=get(handles.zeichnen,'Value');
%save help1 -v6 zeichnen -append;


function help1_Callback(hObject, eventdata, handles)
helpdlg('This routine solves a parameter dependent problem, i.e. the solver provides a solution for each value of the parameter, following a path in the solution/parameter space. Therefore, the routine should be called from an already calculated solution/parameter pair (sol0, p10) (specified in a *.mat file (for example the default value "last.mat")), being a starting point in the path. In each following step the routine modifies both, the solution and the parameter. The value of the parameter p1_0 has to be specified at the very end of the field "Boundary/ Additional conditions", for example in the form p1=7.This routine should only be called with an already calculated solution specified in a *.mat file (for example the default value "last.mat"). The routine modifies the parameter p1 which has to exist in the file! The equation specifying the parameter MUST be in the field "Boundary- / Additional conditions" AFTER the Boundary- / Additional conditions (for example in the form p1=7)','Help');


% --- Executes on button press in help2.
function help2_Callback(hObject, eventdata, handles)
helpdlg('Specify the *.mat file of your initial solution/parameter pair.','Help');


% --- Executes on button press in help2.
function help3_Callback(hObject, eventdata, handles)
helpdlg('Specify the folder where the data, solutions and the corresponding values of parameters should be stored. Provide a name for a file containing a solution/parameter pair, for example solpar.mat. You will then find numbered solpari.mat files containing the name you have chosen. The solpar.mat file (without a number i) contains your path data specified below, see "Pathdata matrix".','Help');

% --- Executes on button press in help4.
function help4_Callback(hObject, eventdata, handles)
helpdlg('In order to follow the path in the solution/parameter space one has to specify the number of steps and the constant stepsize for each step. It may be necessary, when going around turning points, to manually change the stepsize (decrease it). When the run prematurely terminates, because the stepsize was too large, you can return to any successfully computed previous point in the solution/parameter path. Input the number of the solution/parameter pair you want to go back to in "Inital index" and rerun. All following data points will be replaced by your new results.','Help');


% --- Executes on button press in help5.
function help5_Callback(hObject, eventdata, handles)
helpdlg('Initial index of the file-names, cf. previous "?"!','Help');




% --- Executes on button press in help6.
function help6_Callback(hObject, eventdata, handles)
helpdlg('Input an n x 2 matrix, with n being the number of paths you want to show, to save significant data from your pathfollowing strategy. To graphically show the results, one has to input a Matlab matrix in the field "Pathdata matrix". Specify the number of the solution component in the first column (0 represents the parameters). In the second column one has to specify at which point the solution component has to be evaluated (it is not necessary to choose exactly a mesh point, it can be any point in the interval of integration). For example: [3,0.4;0,2;1,0] saves 3 paths, the first one shows z3(0.4) against the pathfollowing parameter p1, the second one shows another parameter denoted by p2 against p1, and the third z1(0) against p1. If one wishes to express the maximum norm of a solution component against p1, then one needs to use a column vector [2;1]. This vec tor means that the maximum norm, 1, of the second solution component, 2, will be plotted against p1.','Help');






function hitpoint_Callback(hObject, eventdata, handles)
% hObject    handle to hitpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hitpoint as text
%        str2double(get(hObject,'String')) returns contents of hitpoint as a double


% --- Executes during object creation, after setting all properties.
function hitpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hitpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help7.
function help7_Callback(hObject, eventdata, handles)
helpdlg('Type a value which you would like to be hit by the path. The algorithm will choose the pathdata in a way that - if passing the hitpoint - one entry will be close to this value.','Help');
