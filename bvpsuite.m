function varargout = bvpsuite(varargin)





% Matlab Code BVPSUITE1.1 (Release July 2010). Authors: Georg Kitzhofer, 
% Othmar Koch, Gernot Pulverer, Christa Simon and Ewa B. Weinmueller. 
% Institute for Analysis and Scientific Computing, Vienna University of Technology, 
% Vienna, Austria. BVPSUITE aims for the efficient numerical solution of boundary 
% value problems in ordinary differential equations. For more information and for 
% the manual contact http://www.math.tuwien.ac.at/~ewa or
% e.weinmueller@tuwien.ac.at.
% Graphical User Interface to BVPSUITE1.1 (Release July 2010)

% Programmed by Georg Kitzhofer 2007

% The Matlab guide-function was used to construct the GUI:

% bvpsuite M-file for bvpsuite.fig
%      bvpsuite, by itself, creates a new bvpsuite or raises the existing
%      singleton*.
%
%      H = bvpsuite returns the handle to a new bvpsuite or the handle to
%      the existing singleton*.
%
%      bvpsuite('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in bvpsuite.M with the given input
%      arguments.
%
%      bvpsuite('Property','Value',...) creates a new bvpsuite or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bvpsuite_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bvpsuite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bvpsuite

% Last Modified by GUIDE v2.5 15-Aug-2009 14:01:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bvpsuite_OpeningFcn, ...
                   'gui_OutputFcn',  @bvpsuite_OutputFcn, ...
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

% --- Executes just before bvpsuite is made visible.
function bvpsuite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bvpsuite (see VARARGIN)

% Choose default command line output for bvpsuite
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%C: set Gaussian
%C: set Initial Value 
%C: set Saved Profile



set(handles.checkbox4,'Value',1);
set(handles.zeichnen,'Value',0);


if isempty(get(handles.edit1,'String'))
    set(handles.radiobutton4,'Value',1);    
    set(handles.automatic,'Value',1);
    set(handles.edit2,'String','1');
    set(handles.radiobutton8,'Value',0);
    set(handles.radiobutton12,'Value',0);
    set(handles.radiobutton1,'Value',1);
    set(handles.checkbox4,'Value',1);
    set(handles.meshadaptation,'Value',1);
    
    if length(strfind(char(cd),'\'))>0
        set(handles.pfad,'String',strcat(char(cd),'\')); %f�r Windows
    else
        set(handles.pfad,'String',strcat(char(cd),'/')); %f�r Linux/Unix
    end


end 


if get(handles.automatic,'Value')==1 %if the number of collocation points should be chosen automatically.
   set(handles.edit13,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   set(handles.edit13,'Enable','off');

load options.mat   
   
%if no bvpfile is open 
if isempty(get(handles.edit1,'String'))
   

  
abstol =1e-10;
reltol = 1e-10;
maxiter= 90000000;
maxfunevals=90000000;
updatejacfactor = 0.5;
lambdamin = 0.001;
switchtoffnfactor = 0.5;
ausgabe =1;
TRM=0;
feinesgitter = 0;
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
   
save options bvpopt abstol reltol maxiter maxfunevals updatejacfactor lambdamin switchtoffnfactor ausgabe TRM abstolgitter reltolgitter K wiederholung feinesgitter n0;
end    
   
 
%Automatic choice of the collocation points depending on the given
%absolute tolerance
  if  min(abstolgitter) <= 10^-6                             

    set(handles.edit13,'String',8)
    
  elseif min(abstolgitter) > 10^-6 && min(abstolgitter) <= 10^-4
      
    set(handles.edit13,'String',6);
  
  elseif min(abstolgitter) > 10^-4 && min(abstolgitter) <= 10^-2
      
    set(handles.edit13,'String',4);     
  elseif min(abstolgitter) > 10^-2 
  
    set(handles.edit13,'String',2);     
  end     
    
    

end







% UIWAIT makes bvpsuite wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bvpsuite_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%varargout{1}=get(handles.edit1,'String');
%varargout{2}=get(handles.pfad,'String');

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%File handling:

function edit2_Callback(hObject, eventdata, handles)
function oeffnen_Callback(hObject, eventdata, handles)

set(handles.checkbox2,'Value',0);
set(handles.radiobutton4,'Value',1);
set(handles.startprofillastmat,'Value',0);

antwort1=get(handles.edit1,'String');
%varargout{1}=antwort1;
%get name for bvpfile!!!!!
bvpfile=char(antwort1);
pfad=char(get(handles.pfad,'String'));
checkm=strfind(bvpfile,'.m');
checkm2=strfind(bvpfile,'last.log');
checkm3=strfind(bvpfile,'bvpsuite.m');



if ((length(checkm)==0) & (length(checkm2)==0))
    err('bvps_errdlg01');
    return;
end
if(length(checkm3)~=0)
    err('bvps_errdlg02');
    return;
end
datei=fopen(strcat(pfad,bvpfile),'r');
if datei==-1
    err('bvps_errdlg03');
    return;
end

%C: l�dt last.log in GUI!! 
if(length(checkm2)~=0) %log file wurde geladen
      if (datei==-1)
          err('bvps_errdlg04');
          return;
      end
       
      inhaltdatei=fread(datei);
      inhaltdatei=char(inhaltdatei');
      fclose(datei);
      help=inhaltdatei(strfind(inhaltdatei,'%dateiname'):strfind(inhaltdatei,'%#dateiname')-1);
      help=strrep(help,'%dateiname','');
      help=help(2:length(help));
      bvpfile=help;
      
       
      help=inhaltdatei(strfind(inhaltdatei,'%Endpoint'):strfind(inhaltdatei,'%#Endpoint')-1);
      help=strrep(help,'%Endpoint','');
      help=help(2:length(help));
      Endpoint=help;
      
      
      
      help=inhaltdatei(strfind(inhaltdatei,'%startwert'):strfind(inhaltdatei,'%#startwert')-1);
      help=strrep(help,'%startwert','');
      help=help(2:length(help));
      startwert=help;
     
      
      help=inhaltdatei(strfind(inhaltdatei,'%dim'):strfind(inhaltdatei,'%#dim')-1);
      help=strrep(help,'%dim','');
      help=help(2:length(help));
      dim=help;
      
 
      % EVP %
      help=inhaltdatei(strfind(inhaltdatei,'%EVP'):strfind(inhaltdatei,'%#EVP')-1);
      help=strrep(help,'%EVP','');
      help=help(2:length(help));
      EVP=help;
      
            
      
      
      % Infinite %
      
      
      help=inhaltdatei(strfind(inhaltdatei,'%Infinite'):strfind(inhaltdatei,'%#Infinite')-1);
      help=strrep(help,'%Infinite','');
      help=help(2:length(help));
      Infinite=help;       
      
      help=inhaltdatei(strfind(inhaltdatei,'%auto'):strfind(inhaltdatei,'%#auto')-1);
      help=strrep(help,'%auto','');
      help=help(2:length(help));
      auto=help;      
      
      
  
            
      help=inhaltdatei(strfind(inhaltdatei,'%standard'):strfind(inhaltdatei,'%#standard')-1);
      help=strrep(help,'%standard','');
      help=help(2:length(help));
      standard=help;
      help=inhaltdatei(strfind(inhaltdatei,'%rho'):strfind(inhaltdatei,'%#rho')-1);
      help=strrep(help,'%rho','');
      help=help(2:length(help));
      rho=help;
      help=inhaltdatei(strfind(inhaltdatei,'%x1'):strfind(inhaltdatei,'%#x1')-1);
      help=strrep(help,'%x1','');
      help=help(2:length(help));
      x1=help;
      help=inhaltdatei(strfind(inhaltdatei,'%g'):strfind(inhaltdatei,'%#g')-1);
      help=strrep(help,'%g','');
      help=help(2:length(help));
      g=help;
      help=inhaltdatei(strfind(inhaltdatei,'%r1'):strfind(inhaltdatei,'%#r1')-1);
      help=strrep(help,'%r1','');
      help=help(2:length(help));
      r1=help;
      help=inhaltdatei(strfind(inhaltdatei,'%parameter'):strfind(inhaltdatei,'%#parameter')-1);
      help=strrep(help,'%parameter','');
      help=help(2:length(help));
      parameter=help;
      
      help=inhaltdatei(strfind(inhaltdatei,'%c'):strfind(inhaltdatei,'%#c')-1);
      help=strrep(help,'%c','');
      help=help(2:length(help));
      c=help;
      
      help=inhaltdatei(strfind(inhaltdatei,'%variablen'):strfind(inhaltdatei,'%#variablen')-1);
      help=strrep(help,'%variablen','');
      help=help(2:length(help));
      variablen=help;
      
      if strfind(inhaltdatei,'%abstol')
          help=inhaltdatei(strfind(inhaltdatei,'%abstol'):strfind(inhaltdatei,'%#abstol')-1);
          help=strrep(help,'%abstol','');
          help=help(3:length(help));
          abstolgitter=str2num(help);
          if length(strfind(char(version),'R14'))==0
              save options -append abstolgitter;
          else
              save options -v6 -append abstolgitter;
          end
      end
     
      if strfind(inhaltdatei,'%reltol')
          help=inhaltdatei(strfind(inhaltdatei,'%reltol'):strfind(inhaltdatei,'%#reltol')-1);
          help=strrep(help,'%reltol','');
          help=help(3:length(help));
          reltolgitter=str2num(help);
          if length(strfind(char(version),'R14'))==0
              save options -append reltolgitter;
          else
              save options -v6 -append reltolgitter;
          end
      end

      
      set(handles.edit1,'String',bvpfile);
      set(handles.pfad,'String',strcat(char(cd),'\'))
      warndlg('The file last.log does not save the directory path!','Warning');
      set(handles.edit2,'String',startwert);
      set(handles.edit11,'String',dim);
      set(handles.edit24,'String',Endpoint);
    % set(handles.automatic,'Value',str2num(auto));
      
        %% changes EVP Christa %%
      
            checkEVP= length(strfind(EVP,'true'));
            
            if checkEVP 
            set(handles.radiobutton8,'Value',1) 
            else
             set(handles.radiobutton8,'Value',0) 
            end 
            
            
            
            
            
            
            
            checkInfinite= length(strfind(Infinite,'true'));
            
            if checkInfinite 
                set(handles.radiobutton12,'Value',1) 
            else
                set(handles.radiobutton12,'Value',0) 
            end 
            
            checkAuto= length(strfind(auto,'true'));
            
            if checkAuto ~= 0
                set(handles.automatic,'Value',1);
                set(handles.edit13,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
                set(handles.edit13,'Enable','off');
            else
                set(handles.automatic,'Value',0);
                set(handles.edit13,'BackgroundColor','white');
                set(handles.edit13,'Enable','on');
            end
            
            
            
            
            
            

            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

      
      checkstand=length(strfind(standard,'true'));
      checkgauss=length(strfind(rho,'[1'));
      checklobatto=length(strfind(rho,'[2'));
      if (checkstand~=0 && checkgauss~=0)  
          set(handles.radiobutton1,'Value',1);
          set(handles.radiobutton2,'Value',0);
          set(handles.radiobutton3,'Value',0);
          set(handles.lobatto,'Value',0);
      end
      if (checkstand~=0 && checklobatto~=0)
          set(handles.radiobutton1,'Value',0);
          set(handles.radiobutton2,'Value',0);
          set(handles.radiobutton3,'Value',0);
          set(handles.lobatto,'Value',1);
      end
      if (checkstand==0)
          set(handles.radiobutton1,'Value',0);
          set(handles.radiobutton2,'Value',1);
          set(handles.radiobutton3,'Value',0);
          set(handles.lobatto,'Value',0);
      end
      
          
          if checkstand ~= 0
              helprho=rho(strfind(rho,' ')+1:strfind(rho,']')-1);
              set(handles.edit13,'String',helprho);
          else
              set(handles.edit13,'String',rho);
          end 
          
          
      
      
      set(handles.edit14,'String',x1);
      set(handles.edit15,'String',g);
      set(handles.edit20,'String',r1);
      set(handles.parameter,'String',parameter);
      set(handles.c,'String',c);
      set(handles.variablen,'String',variablen);
      
      if length(strfind(inhaltdatei,'%lambda'))~=0
          
          help=inhaltdatei(strfind(inhaltdatei,'%lambda'):strfind(inhaltdatei,'%#lambda')-1);
          help=strrep(help,'%lambda','');
          help=help(3:length(help));
          lambda=help;  
          set(handles.lambda,'String',lambda);
      end    
      
      
      if length(strfind(inhaltdatei,'%startprofilgitter'))~=0
          set(handles.checkbox4,'Value',1);
          help=inhaltdatei(strfind(inhaltdatei,'%startprofilgitter'):strfind(inhaltdatei,'%#startprofilgitter')-1);
          help=strrep(help,'%startprofilgitter','');
          help=help(2:length(help));
          startprofilgitter=help;
          
          set(handles.startprofilgitter,'String',startprofilgitter);
          
          help=inhaltdatei(strfind(inhaltdatei,'%startprofilwerte'):strfind(inhaltdatei,'%#startprofilwerte')-1);
          help=strrep(help,'%startprofilwerte','');
          help=help(2:length(help));
          startprofilwerte=help;
          set(handles.startprofilwerte,'String',startprofilwerte);
          
          if length(strfind(inhaltdatei,'%startprofilparameter'))~=0
              help=inhaltdatei(strfind(inhaltdatei,'%startprofilparameter'):strfind(inhaltdatei,'%#startprofilparameter')-1);
              help=strrep(help,'%startprofilparameter','');
              help=help(2:length(help));
              startprofilparameter=help;
              set(handles.startprofilparameter,'String',startprofilparameter);
          else
              set(handles.startprofilparameter,'String','');
          end
      else
          set(handles.checkbox4,'Value',0);
      end

      
  elseif (datei~=-1)
 
      %C:wenn man ein bestehendes file neu l�dt
      
      %C: speichere inhalt von datei in inhaltdatei
      inhaltdatei=fread(datei);
      inhaltdatei=char(inhaltdatei');
     
      
      gueltig=strfind(inhaltdatei,'Kontrollnummer43753976430976655144');
      gueltig2=strfind(inhaltdatei,'Kontrollnummer43753976430976655143');
      gueltig1=strfind(inhaltdatei,'Kontrollnummer43753976430976655145');
      if ((length(gueltig)==0) && (length(gueltig2)==0) && length(gueltig1)==0)
          err('bvps_errdlg05');
          return;
      end
      if (length(gueltig2)>0)
          warndlg('The chosen file is not up to date! Save again! Make sure that if your problem consists of one equation that the solution is denoted by z1 instead of z! Use bvpsuite-help question-marks for further information on special fields!','Warning');
          gueltig1=1;
      elseif (length(gueltig)>0)
          warndlg('The chosen file is not up to date! Save again! Make sure that if your problem consists of one equation that the solution is denoted by z1 instead of z! Use bvpsuite-help question-marks for further information on special fields!','Warning');
          gueltig1=1;
      end
      fclose(datei);
      
      help=inhaltdatei(strfind(inhaltdatei,'%Endpoint'):strfind(inhaltdatei,'%#Endpoint')-1);
      help=strrep(help,'%Endpoint','');
      help=help(3:length(help));
      Endpoint=help;
      
      
      help=inhaltdatei(strfind(inhaltdatei,'%startwert'):strfind(inhaltdatei,'%#startwert')-1);
      help=strrep(help,'%startwert','');
      help=help(3:length(help));
      startwert=help;
      
     
      
      help=inhaltdatei(strfind(inhaltdatei,'%dim'):strfind(inhaltdatei,'%#dim')-1);
      help=strrep(help,'%dim','');
      help=help(3:length(help));
      dim=help;
      
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%% changes EVP Christa %%%%%%%%
     
         
      
      help=inhaltdatei(strfind(inhaltdatei,'%EVP'):strfind(inhaltdatei,'%#EVP')-1);
      help=strrep(help,'%EVP','');
      help=help(3:length(help));
      EVP=help;
      
     
      
      
      help=inhaltdatei(strfind(inhaltdatei,'%Infinite'):strfind(inhaltdatei,'%#Infinite')-1);
      help=strrep(help,'%Infinite','');
     help=help(3:length(help));
      Infinite=help;
      
      
      help=inhaltdatei(strfind(inhaltdatei,'%auto'):strfind(inhaltdatei,'%#auto')-1);
      help=strrep(help,'%auto','');
      help=help(2:length(help));
      auto=help;      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
      help=inhaltdatei(strfind(inhaltdatei,'%standard'):strfind(inhaltdatei,'%#standard')-1);
      help=strrep(help,'%standard','');
      help=help(3:length(help));
      standard=help;
                 
      help=inhaltdatei(strfind(inhaltdatei,'%rho'):strfind(inhaltdatei,'%#rho')-1);
      help=strrep(help,'%rho','');
      help=help(3:length(help));
      rho=help;
      help=inhaltdatei(strfind(inhaltdatei,'%x1'):strfind(inhaltdatei,'%#x1')-1);
      help=strrep(help,'%x1','');
      help=help(3:length(help));
      x1=help;
      help=inhaltdatei(strfind(inhaltdatei,'%g'):strfind(inhaltdatei,'%#g')-1);
      help=strrep(help,'%g','');
      help=help(3:length(help));
      g=help;
      help=inhaltdatei(strfind(inhaltdatei,'%r1'):strfind(inhaltdatei,'%#r1')-1);
      help=strrep(help,'%r1','');
      help=help(3:length(help));
      r1=help;
      help=inhaltdatei(strfind(inhaltdatei,'%parameter'):strfind(inhaltdatei,'%#parameter')-1);
      help=strrep(help,'%parameter','');
      help=help(3:length(help));
      parameter=help;
      help=inhaltdatei(strfind(inhaltdatei,'%c'):strfind(inhaltdatei,'%#c')-1);
      help=strrep(help,'%c','');
      help=help(3:length(help));
      c=help;
      help=inhaltdatei(strfind(inhaltdatei,'%variablen'):strfind(inhaltdatei,'%#variablen')-1);
      help=strrep(help,'%variablen','');
      help=help(3:length(help));
      variablen=help;
      
      if strfind(inhaltdatei,'%abstol')
          help=inhaltdatei(strfind(inhaltdatei,'%abstol'):strfind(inhaltdatei,'%#abstol')-1);
          help=strrep(help,'%abstol','');
          help=help(3:length(help));
          abstolgitter=str2num(help);
          if length(strfind(char(version),'R14'))==0
              save options -append abstolgitter;
          else
              save options -v6 -append abstolgitter;
          end
      end
     
      if strfind(inhaltdatei,'%reltol')
          help=inhaltdatei(strfind(inhaltdatei,'%reltol'):strfind(inhaltdatei,'%#reltol')-1);
          help=strrep(help,'%reltol','');
          help=help(3:length(help));
          reltolgitter=str2num(help);
          if length(strfind(char(version),'R14'))==0
              save options -append reltolgitter;
          else
              save options -v6 -append reltolgitter;
          end
      end
     
       
      %changes EVP
            checkEVP= length(strfind(EVP,'true'));
            
            
            if checkEVP ~= 0
            set(handles.radiobutton8,'Value',1);
            else
             set(handles.radiobutton8,'Value',0); 
            end 

           
            
            checkInfinite= length(strfind(Infinite,'true'));
            
            
            if checkInfinite ~= 0
            set(handles.radiobutton12,'Value',1);
            else
             set(handles.radiobutton12,'Value',0); 
            end 
            
          
            checkAuto= length(strfind(auto,'true'));
            
            
            if checkAuto ~= 0
            set(handles.automatic,'Value',1);
            set(handles.edit13,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
            set(handles.edit13,'Enable','off');
            else
             set(handles.automatic,'Value',0); 
             set(handles.edit13,'BackgroundColor','white');
             set(handles.edit13,'Enable','on');
            end 
            
            
            
      
             set(handles.edit24,'String',Endpoint)
      
   
  
      
      set(handles.edit2,'String',startwert)
      set(handles.edit11,'String',dim);
      checkstand=length(strfind(standard,'true'));
  
      checkgauss=length(strfind(rho,'[1'));
      checklobatto=length(strfind(rho,'[2'));
      if (checkstand~=0 && checkgauss~=0)  
          set(handles.radiobutton1,'Value',1);
          set(handles.radiobutton2,'Value',0);
          set(handles.radiobutton3,'Value',0);
          set(handles.lobatto,'Value',0);
      end
      if (checkstand~=0 && checklobatto~=0)
          set(handles.radiobutton1,'Value',0);
          set(handles.radiobutton2,'Value',0);
          set(handles.radiobutton3,'Value',0);
          set(handles.lobatto,'Value',1);
      end
      if (checkstand==0)
          set(handles.radiobutton1,'Value',0);
          set(handles.radiobutton2,'Value',1);
          set(handles.radiobutton3,'Value',0);
          set(handles.lobatto,'Value',0);
      end
      if checkstand~=0
          helprho=rho(strfind(rho,' ')+1:strfind(rho,']')-1);
          set(handles.edit13,'String',helprho);
      else
          set(handles.edit13,'String',rho);
      end
      set(handles.edit14,'String',x1);
      set(handles.edit15,'String',g);
      set(handles.edit20,'String',r1);
      set(handles.parameter,'String',parameter);
      set(handles.c,'String',c);
      set(handles.variablen,'String',variablen);
      
      
      if length(strfind(inhaltdatei,'%lambda'))~=0
          
           help=inhaltdatei(strfind(inhaltdatei,'%lambda'):strfind(inhaltdatei,'%#lambda')-1);
          help=strrep(help,'%lambda','');
          help=help(3:length(help));
          lambda=help;  
          set(handles.lambda,'String',lambda);
      end    
      
      
      
      
      if length(strfind(inhaltdatei,'%startprofilgitter'))~=0
          set(handles.checkbox4,'Value',1);
          help=inhaltdatei(strfind(inhaltdatei,'%startprofilgitter'):strfind(inhaltdatei,'#startprofilgitter')-1);
          help=strrep(help,'%startprofilgitter','');
          help=help(3:length(help));
          startprofilgitter=help;
          
          set(handles.startprofilgitter,'String',startprofilgitter);
          
          help=inhaltdatei(strfind(inhaltdatei,'%startprofilwerte'):strfind(inhaltdatei,'#startprofilwerte')-1);
          help=strrep(help,'%startprofilwerte','');
          help=help(3:length(help));
          startprofilwerte=help;
          
          set(handles.startprofilwerte,'String',startprofilwerte);
          
          if length(strfind(inhaltdatei,'%startprofilparameter'))~=0
              help=inhaltdatei(strfind(inhaltdatei,'%startprofilparameter'):strfind(inhaltdatei,'%#startprofilparameter')-1);
              help=strrep(help,'%startprofilparameter','');
              help=help(3:length(help));
              startprofilparameter=help;
              set(handles.startprofilparameter,'String',startprofilparameter);
          else
              set(handles.startprofilparameter,'String','');
          end
          
      else
          set(handles.checkbox4,'Value',0);
      end

end

    


function edit11_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit11_Callback(hObject, eventdata, handles)

function edit12_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit12_Callback(hObject, eventdata, handles)


function edit13_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit13_Callback(hObject, eventdata, handles)

function edit14_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit14_Callback(hObject, eventdata, handles)

function edit15_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit15_Callback(hObject, eventdata, handles)
function edit16_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit16_Callback(hObject, eventdata, handles)
function edit17_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit17_Callback(hObject, eventdata, handles)

%********************SAVE*******************************
function speichern_Callback(hObject, eventdata, handles)

try
    test=diff(str2sym('x^2'),'x');
catch
    fprintf(1,'You need the Matlab symbolic toolbox to use this feature!\n');
    fprintf(1,'Check the Matlab help for further information!\n');
    err('bvps_errdlg30');
    error('ERROR');
end



%if get(handles.radiobutton12,'Value') && length(strfind(char(version),'R2009a'))>0
    
%  err('version');
%  return;
    
%end     



antwort1=get(handles.edit1,'String');
bvpfile=char(antwort1);
pfad=char(get(handles.pfad,'String'));

checkm=strfind(bvpfile,'.m');
checkm2=strfind(bvpfile,'last.log');
checkm3=strfind(bvpfile,'bvpsuite.m');
if ((length(checkm)==0) & (length(checkm2)==0))
    err('bvps_errdlg06');
    return;
end
if(length(checkm3)~=0)
    err('bvps_errdlg07');
    return;
end
sicherheit='aaa';
if length(pfad)>0
    datei=fopen(strcat(pfad,bvpfile),'r');
else
    datei=fopen(bvpfile,'r');
end
if datei~=-1
    text=strcat('The file ???',bvpfile,' already exists. Are you sure you want to overwrite it?');
    text=strrep(text,'???','');
    sicherheit=questdlg(text,'Overwrite?');
end
if length(sicherheit)~=3
    return;
end
if datei~=-1
inhaltdatei=fread(datei);
inhaltdatei=char(inhaltdatei');
gueltig=strfind(inhaltdatei,'Kontrollnummer43753976430976655144');
gueltig2=strfind(inhaltdatei,'Kontrollnummer43753976430976655143');
gueltigevp=strfind(inhaltdatei,'Kontrollnummer43753976430976655145');
if ((length(gueltig)==0) && (length(gueltig2)==0) && (length(gueltigevp)==0))
      errordlg('The chosen file is not compatible with bvpsuite!','Error');
      return;
end
if (length(gueltig2)>0)
    warndlg('The file is not up to date. Press save again. Make sure that if your problem consists of one equation that the solution is denoted by z1 instead of z! For further informations see the help messages.','Warning');
    gueltig1=1;
elseif (length(gueltig)>0)
    warndlg('The file is not up to date. Press save again. Make sure that if your problem consists of one equation that the solution is denoted by z1 instead of z!','Warning');
    gueltig1=1;
end
fclose(datei);
end

%name of bvpfile is stored in bvpfile!!!!!!
bvpfile=get(handles.edit1,'String');

antwort1=get(handles.edit2,'String');
antwort2=get(handles.edit11,'String');
if get(handles.radiobutton1,'Value')==1 || get(handles.lobatto,'Value')==1
    antwort3='true';
else
    antwort3='false';
end




if get(handles.radiobutton8,'Value')==1
    antwort_EVP='true';
else
    antwort_EVP='false';
end






if get(handles.radiobutton12,'Value')==1 
    antwort_Infinite='true';

else
    antwort_Infinite='false';
end

if get(handles.automatic,'Value')==1
antwort_auto='true';
else
   antwort_auto='false';
end    


 antwort_endpoint_handle=get(handles.edit24,'String');

 antwort_endpoint= str2num(antwort_endpoint_handle);


if antwort_endpoint < 0
    
    errordlg('Insert a real number >= 0!','Error');
    return;
end 


if get(handles.radiobutton2,'Value')==1
    antwort4=get(handles.edit13,'String');
end
if get(handles.radiobutton1,'Value')==1
        help1=get(handles.edit13,'String');
        
        
        if length(help1)==0
            
            
            if isempty(get(handles.edit13,'String'))
                
            err('bvps_errdlg09');
            return;
            
            end 
        end
        if str2num(help1)>15
            err('bvps_errdlg10');
            return;
        end
        help2=length(strfind(help1,' '));
        if help2>0
            err('bvps_errdlg11');
            return;
        end
        antwort4=strcat('[1 ?',get(handles.edit13,'String'),']');
        antwort4=strrep(antwort4,'?','');
        help=length(strfind(help1,'1'));
        
          
        if (help~=0 && length(help1)==1)
            antwort3='false';
            antwort4='[1/2]';
        end
end

if get(handles.lobatto,'Value')==1
        help1=get(handles.edit13,'String');
        if length(help1)==0
            err('bvps_errdlg09');
            return;
        end
        if str2num(help1)>15
            err('bvps_errdlg12');
            return;
        end
        if str2num(help1)<2
            err('bvps_errdlg13');
            return;
        end
        help2=length(strfind(help1,' '));
        if help2>0
            err('bvps_errdlg11');
            return;
        end
        antwort4=strcat('[2 ?',get(handles.edit13,'String'),']');
        antwort4=strrep(antwort4,'?','');
end


if get(handles.radiobutton3,'Value')==1
    antwort4=strcat('1/',num2str(str2num(get(handles.edit13,'String'))+1),':1/',num2str(str2num(get(handles.edit13,'String'))+1),':',get(handles.edit13,'String'),'/',num2str(str2num(get(handles.edit13,'String'))+1));
    help1=get(handles.edit13,'String');
    if length(help1)==0
        err('bvps_errdlg09');
        return;
    end
    help2=length(strfind(help1,' '));
        if help2>0
            err('bvps_errdlg11');
            return;
        end
        help=length(strfind(help1,'1'));
        if (help~=0 && length(help1)==1)
            antwort4='[1/2]';
        end

end    
    
antwort5=get(handles.edit14,'String');

if isempty(antwort5)
    
    
 errordlg('Insert a mesh for the computations!','Error');   
 return;   
end 

antwort6=get(handles.edit15,'String');

   


if isempty(antwort6)
   errordlg('Fill in the field ''Equations''.!','Error'); 
    return;
end 

%boundary conditions
antwort7=get(handles.edit20,'String');

if isempty(antwort7)
   errordlg('Fill in the field ''Boundary / Additional Conditions''.!','Error'); 
   return;
end 



%�ber anzahl der variablen schleifen


%bvpsuite can only cope with Dirichlet boundary conditions at infinity.
%If the user enters wrong boundary conditions an error message pops up. 
if strcmp(antwort_Infinite,'true')
for m=1:length(antwort2)

        
  l=1;
  while l <= str2num(antwort2(m)) 
     
      
    
    if strfind(antwort7,strcat('z',num2str(m),'''','(b)'))
        
        err('bvps_errdlg15');
        return;
    elseif strfind(antwort7,strcat('z',num2str(m),'d',num2str(l),'(b)'))
        
        err('bvps_errdlg15');
        return;
    
    end
 
    l=l+1;

  end
  

end

end 
    
parameter=get(handles.parameter,'String');
variablen=get(handles.variablen,'String');
c=get(handles.c,'String');
if length(parameter)==0
    parameter='0';
end




if length(antwort2)==0
  
    antwort2= '[]';

end

if length(c)==0
    c='[]';
end
load options
abstol=num2str(abstolgitter);
reltol=num2str(reltolgitter);



datei=fopen('last.log','w');
fprintf(datei,'%%dateiname\n%s%%#dateiname\n%%Endpoint\n%s%%#Endpoint\n%%startwert\n%s%%#startwert\n',bvpfile,num2str(antwort_endpoint),char(antwort1));
fprintf(datei,'%%dim\n%s%%#dim\n%%parameter\n%s%%#parameter\n%%c\n%s%%#c\n%%variablen\n%s%%#variablen\n%%standard\n%s%%#standard\n%%EVP\n%s%%Infinite\n%s%%#Infinite\n',char(antwort2),char(parameter),char(c),char(variablen),char(antwort3),char(antwort_EVP),char(antwort_Infinite));  

fprintf(datei,'%%auto\n%s%%#auto\n',char(antwort_auto));

fprintf(datei,'%%rho\n%s%%#rho\n%%x1\n%s%%#x1\n%%abstol\n%s%%#abstol\n%%reltol\n%s%%#reltol\n',char(antwort4),char(antwort5),abstol,reltol);
antwort6=strrep(strrep(antwort6,'[',''),']','');
antwort7=strrep(strrep(antwort7,'[',''),']','');



fprintf(datei,'%%g\n%s%%#g\n%%r1\n%s%%#r1\n',char(antwort6),char(antwort7));

if get(handles.checkbox4,'Value')==1
    fprintf(datei,'%%startprofilgitter\n%s%%#startprofilgitter\n',get(handles.startprofilgitter,'String'));
    fprintf(datei,'%%startprofilwerte\n%s%%#startprofilwerte\n',get(handles.startprofilwerte,'String'));
    
end
if length(str2num(get(handles.startprofilparameter,'String')))>0
    fprintf(datei,'%%startprofilparameter\n%s%%#startprofilparameter\n',get(handles.startprofilparameter,'String'));
end
if length(str2num(get(handles.lambda,'String')))>0
    fprintf(datei,'%%lambda\n%s%%#lambda\n',get(handles.lambda,'String'));
end

fclose(datei);

antwort6=strcat('[',char(antwort6),']');
antwort7=strcat('[',char(antwort7),']');

antwort2_old=antwort2;
antwort6_old=antwort6;
antwort7_old=antwort7;

%Ersetzen der Variablen in den Gleichungen (antwort6)
variablenx=regexprep(variablen,'\s*','');
if length(variablenx)>0
    if variablenx(length(variablenx))~=';'
        variablenx=strcat(variablenx,';');
    end
end
antwort6x=antwort6;
antwort7x=antwort7;
while(length(strfind(variablenx,'='))>0)
    gleichheitszeichenstelle=strfind(variablenx,'=');
    strichpunktstelle=strfind(variablenx,';');
    variablenname=variablenx(1:gleichheitszeichenstelle(1)-1);
    variablenwert=variablenx(gleichheitszeichenstelle(1)+1:strichpunktstelle(1)-1);
    antwort6x = regexprep(antwort6x,strcat('(?<=^|[^a-zA-Z0-9_])',variablenname,'(?=$|[^a-zA-Z0-9_])'),strcat('(',variablenwert,')'));
    antwort7x = regexprep(antwort7x,strcat('(?<=^|[^a-zA-Z0-9_])',variablenname,'(?=$|[^a-zA-Z0-9_])'),strcat('(',variablenwert,')'));
    variablenx = regexprep(variablenx(strichpunktstelle(1)+1:length(variablenx)),strcat('(?<=^|[^a-zA-Z0-9_])',variablenname,'(?=$|[^a-zA-Z0-9_])'),strcat('(',variablenwert,')'));
end






%******************************************************
%      Call EVP module!!!!
%******************************************************



 if strcmp(antwort_EVP,'true')
  
   [antwort2,antwort6x,antwort7x] = EVPmoduleEdit(antwort2,antwort6x,antwort7x);    
     
   
 end

 
%******************************************************
%      Call Transformation module!!!!
%******************************************************

if strcmp(antwort_Infinite,'true')
    
   
     anz_glei=length(str2num(char(antwort2)));
        
    [antwort2,antwort6x,antwort7x] = trafomodule(antwort6x,antwort7x,anz_glei,antwort_EVP,antwort2,antwort_endpoint,parameter); 
    
    
end

mfileschreiben(antwort1,antwort2,antwort3,antwort4,antwort5,antwort6x,antwort7x,antwort_EVP,antwort_Infinite,antwort_endpoint,bvpfile,pfad,parameter,c,antwort_auto)
if length(pfad)>0
    datei=fopen(strcat(pfad,bvpfile),'a');
else
    datei=fopen(bvpfile,'a');
end

if strcmp(antwort_EVP,'true') || strcmp(antwort_Infinite,'true')
    
    antwort2=antwort2_old ;
    antwort6 = antwort6_old;
    antwort7= antwort7_old;
   
end




fprintf(datei,'%%Values read by bvpsuite GUI:\n');
fprintf(datei,'%%Endpoint\n%%%s%%#Endpoint\n%%startwert\n%%%s%%#startwert\n',num2str(antwort_endpoint),char(antwort1));
fprintf(datei,'%%dim\n%%%s%%#dim\n%%parameter\n%%%s%%#parameter\n%%c\n%%%s%%#c\n%%variablen\n%%%s%%#variablen\n%%standard\n%%%s%%#standard\n%%Infinite\n%%%s%%#Infinite\n%%EVP\n%%%s%%#EVP\n',char(antwort2),char(parameter),char(c),char(variablen),char(antwort3),char(antwort_Infinite),char(antwort_EVP));
fprintf(datei,'%%rho\n%%%s%%#rho\n%%x1\n%%%s%%#x1\n%%abstol\n%%%s%%#abstol\n%%reltol\n%%%s%%#reltol\n',char(antwort4),char(antwort5),abstol,reltol);
fprintf(datei,'%%auto\n%%%s%%#auto\n',antwort_auto);
antwort6=strrep(strrep(antwort6,'[',''),']','');
antwort7=strrep(strrep(antwort7,'[',''),']','');
fprintf(datei,'%%g\n%%%s%%#g\n%%r1\n%%%s%%#r1\n',char(antwort6),char(antwort7));
if get(handles.checkbox4,'Value')==1
    fprintf(datei,'%%startprofilgitter\n%%%s#startprofilgitter\n',get(handles.startprofilgitter,'String'));
    fprintf(datei,'%%startprofilwerte\n%%%s#startprofilwerte\n',get(handles.startprofilwerte,'String'));
end
if length(str2num(get(handles.startprofilparameter,'String')))>0
    fprintf(datei,'%%startprofilparameter\n%%%s%%#startprofilparameter\n',get(handles.startprofilparameter,'String'));
end
if length(str2num(get(handles.lambda,'String')))>0
    fprintf(datei,'%%lambda\n%%%s%%#lambda\n',get(handles.lambda,'String'));
end
fclose(datei);

fprintf(1,'The file %s has been written!\n',bvpfile);
fprintf(1,'Check the inputs in the file!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deleting function cache (i think this is unnecessary M.Z.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rehash;

helpdlg('The file has been written!','SUCCESS');





function ret=mfileschreiben(antwort1,antwort2,antwort3,antwort4,antwort5,antwort6,antwort7,antwort_EVP,antwort_Infinite,antwort_endpoint,bvpfile,pfad,parameter,c,antwort_auto);
dimension=num2str(length(str2num(char(antwort2))));
ordnung=str2num(char(antwort2));
parameter=str2num(char(parameter));
c=str2num(char(c));
g=char(antwort6);
r1=char(antwort7);


 if strcmp(antwort_Infinite,'true')
 
  
% antwort5 = str2num(antwort5)   
 if antwort_endpoint ~= 0
    trafo_endpoint=1/antwort_endpoint;
    %help=strcat('linspace(0,',num2str(trafo_endpoint)) ;
    %antwort5 = strcat('linspace(0,1,',antwort5,')')
     help=strcat('linspace(0,',num2str(1)) ;   
 else
    help=strcat('linspace(0,',num2str(1)) ;    
     
 end 
 antwort5 = strcat(help,',',antwort5,')');
 
 end     
    
if length(pfad)>0
    datei=fopen(strcat(pfad,bvpfile),'w');
else
    datei=fopen(bvpfile,'w');
end
help=strrep(bvpfile,'.m','');




%C: first line in bvpfile is function definition:
fprintf(datei,'function [ret]=%s(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)\n\n',help);
fprintf(datei,'%%Inputs for bvpsuite\n%%Programmed by Georg Kitzhofer, July 2004\n\n');
fprintf(datei,'%%Kontrollnummer43753976430976655145\n');
fprintf(datei,'if (nargin<15) zc=[]; if (nargin<14) DpRpar=[]; if (nargin<13) DpRgl=[]; if (nargin<12) Dppar=[]; if (nargin<11) Dpgl=[]; if (nargin<10) p=[]; if (nargin<9) D2=[]; if (nargin<8) D1=[]; if (nargin<7) zb=[]; if (nargin<6) za=[]; if (nargin<5) z=[]; if (nargin<4) t=[]; if (nargin<3) R=[];\nif (nargin<2) D=[];end;end;end;end;end;end;end;end;end;end;end;end;end;end;');
fprintf(datei,'\nret=intern(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc);\n\nfunction ret=intern(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)\n');
fprintf(datei,'switch nr\n     case ''x0''\n         for j=1:(length(intern(''x1''))-1)*(intern(''n'')*intern(''m'')+sum(intern(''ordnung'')))+intern(''parameter'')\n');

fprintf(datei,'             u(j,1)=%s;\n',char(antwort1));

fprintf(datei,'         end\n         ret=u;\n     case ''ordnung''\n         ret=[%s];\n     case ''parameter''\n         ret=%s;\n     case ''c''\n         ret=[%s];\n     case ''n''\n         ret=%s;\n',num2str(ordnung),num2str(parameter),num2str(c),char(dimension));
    
fprintf(datei,'     case ''Infinite''\n         ret=%s;\n',char(antwort_Infinite));    
fprintf(datei,'     case ''EVP''\n         ret=%s;\n',char(antwort_EVP));
if strcmp(antwort_Infinite,'true')
 

fprintf(datei,'     case ''Endpoint''\n         ret=%s;\n',char(num2str(antwort_endpoint)));

else
      
  fprintf(datei,'     case ''Endpoint''\n         ret=[];\n');  
    
end 
fprintf(datei,'     case ''standard''\n         ret=%s;\n',char(antwort3));
fprintf(datei,'     case ''rho''\n         ret=%s;\n',char(antwort4));
fprintf(datei,'     case ''x1''\n         ret=%s;\n',char(antwort5));
fprintf(datei,'     case ''g''\n');
gausgabe=prepausgabe(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
fprintf(datei,'         ret=%s;\n',gausgabe);
fprintf(datei,'     case ''Dg''\n         switch D\n');
gsymfunktion=prepdiff(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,ordnung,parameter);
linear=1;
for oi=0:max(ordnung)
    fprintf(datei,'          case %s\n',int2str(oi+1));
  %  helpy=('');
    fprintf(datei,'                 ret=[');
    for ni1=1:str2num(char(dimension))
        for nicounter=1:str2num(char(dimension))
            respect=strcat('z',int2str(nicounter),'d',int2str(oi));
            helpx=char(diff(gsymfunktion(ni1),respect));
	        helpx = regexprep(helpx, ' ', '');
            helpx=prepausgabegdiff(helpx,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
            fprintf(datei,'%s ',helpx);
            if length(strfind(helpx,'z')) || length(strfind(helpx,'p'))>0
                linear=0;
            end
        end
        fprintf(datei,';');
    end
    fprintf(datei,'];\n');
end
fprintf(datei,'         end\n');
fprintf(datei,'     case ''Dpg''\n');
%fprintf(datei,'         switch Dpgl\n');
if length(ordnung)>0
    fprintf(datei,'         ret=[');
end
for ni=1:length(ordnung)
    for pii=1:parameter
        respect=strcat('p',int2str(pii));
        help=char(diff(gsymfunktion(ni),respect));
        help = regexprep(help, ' ', '');
        help=prepausgabegdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
        fprintf(datei,'%s ',help);
        if length(strfind(help,'z')) || length(strfind(help,'p'))>0
            linear=0;
        end
    end
    if parameter>0
        fprintf(datei,';');
    else
        fprintf(datei,'0;');
    end
end
if length(ordnung)>0
    fprintf(datei,'];');
else
    fprintf(datei,'         ret=[];');
end
fprintf(datei,'\n');


fprintf(datei,'                 %%Additional conditions\n');
if length(c)==0
    fprintf(datei,'                 %%za(Komponente,Ableitung)\n');
else
    fprintf(datei,'                 %%zc(Komponente,Ableitung,Intervallstelle c_i)\n');
end
fprintf(datei,'     case ''R''\n');%         switch R\n');
if length(c)==0  
    r1ausgabe=prepausgaber(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
   else
    
    r1ausgabe=prepausgaber_c(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
end
if length(r1ausgabe)==2
    if min(r1ausgabe=='[]')
        r1ausgabe='[]';
    end
end
fprintf(datei,'         ret=%s;\n',r1ausgabe);
fprintf(datei,'     case ''DR''\n          switch R\n');
if length(c)==0
    r1symfunktion=prepdiffr(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
   
    
    
    %Jede Randbedingung hat ihre Ableitung
    
    
    for p=1:sum(ordnung)+parameter

        fprintf(datei,'             case %i\n',p);
        fprintf(datei,'                 switch D\n');

        casew=0;
        for oi=0:max(ordnung)-1
            casew=casew+1;fprintf(datei,'                      case ''a%s''\n',int2str(oi+1));
           fprintf(datei,'                          ret=[');
            %Die �u�ere Schleife bestimmt, welche Funktionen abgeleitet werden
            for nicounter=1:str2num(char(dimension))
                
                
                respect=strcat('za',int2str(nicounter),'d',int2str(oi));
                help=char(diff(r1symfunktion(p),respect));
                help = regexprep(help, ' ', '');
                help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
                fprintf(datei,'%s ',help);
                
                if length(strfind(help,'z'))>0 || length(strfind(help,'p'))>0
                    linear=0;
                end
            end
            fprintf(datei,';');
            fprintf(datei,'];\n',int2str(oi+1));
        end
        casew=0;
        for oi=0:max(ordnung)-1
            casew=casew+1;fprintf(datei,'                      case ''b%s''\n',int2str(oi+1));
            fprintf(datei,'                          ret=[');
            for nicounter=1:str2num(char(dimension))
                respect=strcat('zb',int2str(nicounter),'d',int2str(oi));
                help=char(diff(r1symfunktion(p),respect));
		        help = regexprep(help, ' ', '');
                help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
                fprintf(datei,'%s ',help);
                if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                   linear=0;
                end
            end
            fprintf(datei,';');
            fprintf(datei,'];\n',int2str(oi+1));
        end

        fprintf(datei,'                  end\n');
    end
    fprintf(datei,'          end\n');
    fprintf(datei,'     case ''DpR''\n');
    %fprintf(datei,'         switch DpRgl\n');
    if sum(ordnung)+parameter>0
        fprintf(datei,'         ret=[');
    end
    for ni=1:sum(ordnung)+parameter
        for pii=1:parameter
            respect=strcat('p',int2str(pii));
            help=char(diff(r1symfunktion(ni),respect));
            help = regexprep(help, ' ', '');
            help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
            fprintf(datei,'%s ',help);
            if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                linear=0;
            end
        end
        if parameter>0
            fprintf(datei,';');
        else
            fprintf(datei,'0;');
        end
    end
    if sum(ordnung)+parameter>0
        fprintf(datei,'];');
    else
        fprintf(datei,'         ret=[];');
    end
    fprintf(datei,'\n');
else
    r1symfunktion=prepdiffr_c(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
    %Jede Randbedingung hat ihre Ableitung
    for p=1:sum(ordnung)+parameter

        fprintf(datei,'        case %i\n',p);
        fprintf(datei,'            switch D\n');
        for ci=1:length(c)
            casew=0;
            for oi=0:max(ordnung)-1
                casew=casew+1;fprintf(datei,'                      case ''c%s_%s''\n',int2str(ci),int2str(oi+1));
                fprintf(datei,'                          ret=[');
                %Die �u�ere Schleife bestimmt, welche Funktionen abgeleitet werden
                for nicounter=1:str2num(char(dimension))
                    respect=strcat('zc',num2str(ci),'_',int2str(nicounter),'d',int2str(oi));
                    help=char(diff(r1symfunktion(p),respect));
                    help = regexprep(help, ' ', '');
                    help=prepausgaberdiff_c(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
                    fprintf(datei,'%s ',help);
                    if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                        linear=0;
                    end
                end
                fprintf(datei,';');
                fprintf(datei,'];\n',int2str(oi+1));
            end
        end
        
        
        fprintf(datei,'                  end\n');
    end
    fprintf(datei,'          end\n');
    fprintf(datei,'     case ''DpR''\n');
    %fprintf(datei,'         switch DpRgl\n');
    if sum(ordnung)+parameter>0
        fprintf(datei,'         ret=[');
    end
    for ni=1:sum(ordnung)+parameter
        for pii=1:parameter
            respect=strcat('p',int2str(pii));
            help=char(diff(r1symfunktion(ni),respect));
            help = regexprep(help, ' ', '');
            if length(c)==0
                help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
                if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                    linear=0;
                end
            else
                help=prepausgaberdiff_c(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
                if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                    linear=0;
                end
            end
            fprintf(datei,'%s ',help);
        end
        if parameter>0
            fprintf(datei,';');
        else
            fprintf(datei,'0;');
        end
    end
    if sum(ordnung)+parameter>0
        fprintf(datei,'];');
    else
        fprintf(datei,'         ret=[];');
    end
    fprintf(datei,'\n');
end


fprintf(datei,'     case ''m''\n         if (intern(''standard''))\n             help=intern(''rho'');\n');
fprintf(datei,'             ret=help(2);\n         else\n             ret=length(intern(''rho''));\n');
fprintf(datei,'         end\n     case ''linear''\n         ret=');
if linear
    fprintf(datei,'1;\n end\n');
else
    fprintf(datei,'0;\n end\n');
end
fclose(datei);

%C: ende der funktion mfileschreiben



function ret=prepdiff(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,ordnung,parameter)

%Erkennen der vorkommenden z und markieren durch den String xyq...
for oi=max(ordnung):-1:0
    for ni=length(ordnung):-1:1
        switch oi
            case 0
                help=strcat('z',int2str(ni));
            case 1
                help=strcat('z',int2str(ni),'''');
            case 2
                help=strcat('z',int2str(ni),'''''');
            
        end
        helpneu=strcat('z',int2str(ni),'d',int2str(oi));
        strcat('xyq(',int2str(ni),',',int2str(oi+1),')');
        g=strrep(g,helpneu,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
       
        
         if oi<=2
            g=strrep(g,help,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
        end
    end
end


for oi=max(ordnung):-1:0
    switch oi
        case 2
            help='z''''';
        case 1
            help='z''';
        case 0
            help='z';
       
     
    end
    helpneu=strcat('zd',int2str(oi));
    help2=strcat('xyq(1,',int2str(oi+1),')');
    g=strrep(g,helpneu,help2);
    
    
    if oi<=2
        g=strrep(g,help,help2);
    end
end
%Erkennen der unbekannten Parameter p und ersetzen durch den String yxq
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    g=strrep(g,help,strcat('yxq(',int2str(pii),')'));
end
g=strrep(g,'xyq','z');
g=strrep(g,'yxq','p');
g=strrep(g,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(g,'=');
gbak=g;
for j=1:length(gleichheitszeichen)
    rest=gbak(gleichheitszeichen(j):length(gbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        g=strcat(g(1:stelle-1),')',g(stelle:length(g)));
    else
        g=strrep(g,']',')]');
    end
end
g=strrep(g,'=','-(');
%Ende Gleichheitszeichen

g=eval(strcat('inline(''',g,''',''z'',''p'',''t'')'));

z=sym(0);
for oi=0:max(ordnung)
    for ni=1:length(ordnung)
        help=strcat('z',int2str(ni),'d',int2str(oi));
        z(ni,oi+1)=str2sym(help);
    end
end
p=sym(0);
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    p(pii)=str2sym(help);
end
try
    g=g(z,p,str2sym('t'));
catch
    err('bvps_errdlg31');
    g=g(z,p,str2sym('t'));
end
ret=g;

function ret=prepdiffr(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)

for oi=o-1:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                helpa=strcat('z',int2str(ni),'(a)');
                helpb=strcat('z',int2str(ni),'(b)');
            case 1
                helpa=strcat('z',int2str(ni),'''(a)');
                helpb=strcat('z',int2str(ni),'''(b)');
           case 2
               helpa=strcat('z',int2str(ni),'''''(a)');
               helpb=strcat('z',int2str(ni),'''''(b)');
                
        
        end
        helpaneu=strcat('z',int2str(ni),'d',int2str(oi),'(a)');
        helpbneu=strcat('z',int2str(ni),'d',int2str(oi),'(b)');
        r=strrep(r,helpaneu,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
        r=strrep(r,helpbneu,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
      
        if oi<=2
            r=strrep(r,helpa,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
            r=strrep(r,helpb,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
        end
    end
end
for oi=o-1:-1:0
    
    
    switch oi
        case 1
            helpa='z''(a)';
            helpb='z''(b)';
        case 0
            helpa='z(a)';
            helpb='z(b)';
        case 2
            helpa='z''''(a)';
            helpb='z''''(b)';
    
    end
    helpaneu=strcat('zd',int2str(oi),'(a)');
    helpbneu=strcat('zd',int2str(oi),'(b)');
    help2a=strcat('xyqa(1,',int2str(oi+1),')');
    help2b=strcat('xyqb(1,',int2str(oi+1),')');
    r=strrep(r,helpaneu,help2a);
    r=strrep(r,helpbneu,help2b);
  
    if oi<=2
        r=strrep(r,helpa,help2a);
        r=strrep(r,helpb,help2b);
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end

r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(r,'=');
rbak=r;
for j=1:length(gleichheitszeichen)
    rest=rbak(gleichheitszeichen(j):length(rbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        r=strcat(r(1:stelle-1),')',r(stelle:length(r)));
    else
        r=strrep(r,']',')]');
    end
end
r=strrep(r,'=','-(');
%Ende Gleichheitszeichen



r=eval(strcat('inline(''',r,''',''za'',''zb'',''p'')'));

za=sym(0);
zb=sym(0);
for oi=0:o
    for ni=1:str2num(char(dimension))
        helpa=strcat('za',int2str(ni),'d',int2str(oi));
        helpb=strcat('zb',int2str(ni),'d',int2str(oi));
        za(ni,oi+1)=str2sym(helpa);
        zb(ni,oi+1)=str2sym(helpb);
    end
end
p=sym(0);
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    p(pii)=str2sym(help);
end

try
    r=r(za,zb,p);
catch
    err('bvps_errdlg32');
    r=r(za,zb,p);
end
ret=r;



function ret=prepdiffr_c(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter,c)

for oi=o-1:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                for ci=1:length(c)
                    help2=strcat('z',int2str(ni),'(c',num2str(ci),')');
                    help(ci,1:length(help2))=help2;
                    helplen(ci)=length(help2);
                end
            case 1
                for ci=1:length(c)
                    help2=strcat('z',int2str(ni),'''(c',num2str(ci),')');
                    help(ci,1:length(help2))=help2;
                    helplen(ci)=length(help2);
                end
        end
        for ci=1:length(c)
            help2=strcat('z',int2str(ni),'d',int2str(oi),'(c',num2str(ci),')');
            helpneu(ci,1:length(help2))=help2;
            helpneulen(ci)=length(help2);
        end
        for ci=length(c):-1:1
            r=strrep(r,helpneu(ci,1:helpneulen(ci)),strcat('xyqc(',int2str(ni),',',int2str(oi+1),',',int2str(ci),')'));
        end    
        if oi<=1
            for ci=length(c):-1:1
                r=strrep(r,help(ci,1:helplen(ci)),strcat('xyqc(',int2str(ni),',',int2str(oi+1),',',int2str(ci),')'));
            end
        end
    end
end
for oi=o-1:-1:0
    switch oi
        case 1
            for ci=1:length(c)
                help2=strcat('z''(c',num2str(ci),')');
                help(ci,1:length(help2))=help2;
                helplen(ci)=length(help2);
            end
        case 0
            for ci=1:length(c)
                help2=strcat('z(c',num2str(ci),')');
                help(ci,1:length(help2))=help2;
                helplen(ci)=length(help2);
            end       
    end
    for ci=1:length(c)
        help2=strcat('zd',int2str(oi),'(c',num2str(ci),')');
        helpneu(ci,1:length(help2))=help2;
        helpneulen(ci)=length(help2);
    end
    for ci=1:length(c)
        help2=strcat('xyqc(1,',int2str(oi+1),',',int2str(ci),')');
        help3(ci,1:length(help2))=help2;
        help3len(ci)=length(help2);
    end
    for ci=length(c):-1:1
        r=strrep(r,helpneu(ci,1:helpneulen(ci)),help3(ci,1:help3len(ci)));
    end
    if oi<=1
        for ci=length(c):-1:1
            r=strrep(r,help(ci,1:helplen(ci)),help3(ci,1:help3len(ci)));
        end
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end

r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(r,'=');
rbak=r;
for j=1:length(gleichheitszeichen)
    rest=rbak(gleichheitszeichen(j):length(rbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        r=strcat(r(1:stelle-1),')',r(stelle:length(r)));
    else
        r=strrep(r,']',')]');
    end
end
r=strrep(r,'=','-(');
%Ende Gleichheitszeichen



help2=strcat('inline(''',r,''',''zc'',''p'')');



r=eval(help2)

%Da es nicht m�glich ist, einen Namen zu generieren, der als Variablenname
%weiterverwendet werden kann, wird hier ein zuf�lliges File erzeugt, das
%dann alle notwendigen Eintr�ge in Stringform enth�lt, aufgerufen und
%gleich wieder gel�scht wird.



zc=sym(0);

for oi=0:o
    for ni=1:str2num(char(dimension))
        for ci=1:length(c)
             help2=strcat('zc',num2str(ci),'_',int2str(ni),'d',int2str(oi));
             helpc(ci,1:length(help2))=help2;
             helpclen(ci)=length(help2);
        end
   %     helpb=strcat('zb',int2str(ni),'d',int2str(oi));
        for ci=1:length(c)
            zc(ni,oi+1,ci)=str2sym(helpc(ci,1:helpclen(ci)));
        end
   %     zb(ni,oi+1)=sym(helpb);
    end;
end;
p=sym(0);
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    p(pii)=str2sym(help);
end;
try
    r=r(zc,p);
catch
    err('bvps_errdlg32');
    r=r(zc,p);
end
ret=r;

function ret=prepausgabe(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)

for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                help=strcat('z',int2str(ni));
            case 1
                help=strcat('z',int2str(ni),'''');
            case 2
                help=strcat('z',int2str(ni),'''''');
           
        end
        helpneu=strcat('z',int2str(ni),'d',int2str(oi));
       g=strrep(g,helpneu,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
       if oi<=2      
        g=strrep(g,help,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
       end
    end
end
for oi=o:-1:0
    switch oi       
        case 2
            help='z''''';
        case 1
            help='z''';
        case 0
            help='z';
    end
    helpneu=strcat('zd',int2str(oi));
    help2=strcat('xyq(1,',int2str(oi+1),')');
    g=strrep(g,helpneu,help2);
   
   
    if oi<=2
        g=strrep(g,help,help2);
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    g=strrep(g,help,strcat('yxq(',int2str(pii),')'));
end

g=strrep(g,'xyq','z');
g=strrep(g,'yxq','p');
g=strrep(g,' ','');

%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(g,'=');
gbak=g;
for j=1:length(gleichheitszeichen)
    rest=gbak(gleichheitszeichen(j):length(gbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        g=strcat(g(1:stelle-1),')',g(stelle:length(g)));
    else
        g=strrep(g,']',')]');
    end
end
g=strrep(g,'=','-(');
%Ende Gleichheitszeichen


ret=g;

% Ende prepausgabe 

function ret=prepausgaber(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)

for oi=o-1:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                helpa=strcat('z',int2str(ni),'(a)');
                helpb=strcat('z',int2str(ni),'(b)');
            case 1
                helpa=strcat('z',int2str(ni),'''(a)');
                helpb=strcat('z',int2str(ni),'''(b)');
            case 2
                helpa=strcat('z',int2str(ni),'''''(a)');
                helpb=strcat('z',int2str(ni),'''''(b)');     
                
                
                
        end
        helpaneu=strcat('z',int2str(ni),'d',int2str(oi),'(a)');
        helpbneu=strcat('z',int2str(ni),'d',int2str(oi),'(b)');
        r=strrep(r,helpaneu,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
        r=strrep(r,helpbneu,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
        
        if oi<=2
            r=strrep(r,helpa,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
            r=strrep(r,helpb,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
        end
    end
end
for oi=o-1:-1:0
    switch oi
        case 1
            helpa='z''(a)';
            helpb='z''(b)';
        case 0
            helpa='z(a)';
            helpb='z(b)';
        case 2 
            helpa='z''''(a)';
            helpb='z''''(b)';
            
            
    end
    helpaneu=strcat('zd',int2str(oi),'(a)');
    helpbneu=strcat('zd',int2str(oi),'(b)');
    r=strrep(r,helpaneu,strcat('xyqa(1,',int2str(oi+1),')'));
    r=strrep(r,helpbneu,strcat('xyqb(1,',int2str(oi+1),')'));
    
    if oi<=2
        r=strrep(r,helpa,strcat('xyqa(1,',int2str(oi+1),')'));
        r=strrep(r,helpb,strcat('xyqb(1,',int2str(oi+1),')'));
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end
r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(r,'=');
rbak=r;
for j=1:length(gleichheitszeichen)
    rest=rbak(gleichheitszeichen(j):length(rbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        r=strcat(r(1:stelle-1),')',r(stelle:length(r)));
    else
        r=strrep(r,']',')]');
    end
end
r=strrep(r,'=','-(');
%Ende Gleichheitszeichen
    
    
ret=r;

function ret=prepausgaber_c(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter,c)

for oi=o-1:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                for ci=1:length(c)
                    help2=strcat('z',int2str(ni),'(c',int2str(ci),')');
                    help(ci,1:length(help2))=help2;
                    helplen(ci)=length(help2);
                end
            case 1
                for ci=1:length(c)
                    help2=strcat('z',int2str(ni),'''(c',int2str(ci),')');
                    help(ci,1:length(help2))=help2;
                    helplen(ci)=length(help2);
                end
        end
        for ci=1:length(c)
            help2=strcat('z',int2str(ni),'d',int2str(oi),'(c',int2str(ci),')');
            helpneu(ci,1:length(help2))=help2;
            helpneulen(ci)=length(help2);
        end
        for ci=length(c):-1:1
            r=strrep(r,helpneu(ci,1:helpneulen(ci)),strcat('xyqc(',int2str(ni),',',int2str(oi+1),',',int2str(ci),')'));
        end
        if oi<=1
            for ci=length(c):-1:1
                r=strrep(r,help(ci,1:helplen(ci)),strcat('xyqc(',int2str(ni),',',int2str(oi+1),',',int2str(ci),')'));
            end
        end
    end
end
for oi=o-1:-1:0
    switch oi
        case 1
            for ci=1:length(c)
                help2=strcat('z''(c',int2str(ci),')');
                help(ci,1:length(help2))=help2;
                helplen(ci)=length(help2);
            end
        case 0
            for ci=1:length(c)
                help2=strcat('z(c',int2str(ci),')');
                help(ci,1:length(help2))=help2;
                helplen(ci)=length(help2);
            end
    end
    for ci=length(c):-1:1
        help2=strcat('zd',int2str(oi),'(c',int2str(ci),')');
        helpneu(ci,1:length(help2))=help2;
        helpneulen(ci)=length(help2);
    end
    for ci=length(c):-1:1
        r=strrep(r,helpneu(ci,1:helpneulen(ci)),strcat('xyqc(1,',int2str(oi+1),',',int2str(ci),')'));
    end
    if oi<=1
        for ci=length(c):-1:1
            r=strrep(r,help(ci,1:helplen(ci)),strcat('xyqc(1,',int2str(oi+1),',',int2str(ci),')'));
        end
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end
r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(r,'=');
rbak=r;
for j=1:length(gleichheitszeichen)
    rest=rbak(gleichheitszeichen(j):length(rbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        r=strcat(r(1:stelle-1),')',r(stelle:length(r)));
    else
        r=strrep(r,']',')]');
    end
end
r=strrep(r,'=','-(');
%Ende Gleichheitszeichen
    
    
ret=r;



function ret=prepausgabegdiff(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)
for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        help=strcat('z',int2str(ni),'d',int2str(oi));
        g=strrep(g,help,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    g=strrep(g,help,strcat('yxq(',int2str(pii),')'));
end
g=strrep(g,'xyq','z');
g=strrep(g,'yxq','p');
g=strrep(g,'abs(','abs1(');
ret=g;

function ret=prepausgaberdiff(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)
for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        helpa=strcat('za',int2str(ni),'d',int2str(oi));
        helpb=strcat('zb',int2str(ni),'d',int2str(oi));
        r=strrep(r,helpa,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
        r=strrep(r,helpb,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end
r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,'abs(','abs1(');
ret=r;

function ret=prepausgaberdiff_c(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter,c)
for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        for ci=1:length(c)
            help2=strcat('zc',int2str(ci),'_',int2str(ni),'d',int2str(oi));
            help(ci,1:length(help2))=help2;
            helplen(ci)=length(help2);
        end
        for ci=length(c):-1:1
            r=strrep(r,help(ci,1:helplen(ci)),strcat('xyqc(',int2str(ni),',',int2str(oi+1),',',int2str(ci),')'));
        end
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end
r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,'abs(','abs1(');
ret=r;



%Ende bvpfile generator






function pushbutton7_Callback(hObject, eventdata, handles)

aktuellesverzeichnis=cd;
arbeitsverzeichnis=get(handles.pfad,'String');
cd(arbeitsverzeichnis);

[FileName,PathName,index] = uigetfile('*.m;*.log','Choose the m-file');

if index>0
    set(handles.edit1,'String',FileName);
    set(handles.pfad,'String',PathName);
end
cd(aktuellesverzeichnis);
if index>0
    oeffnen_Callback(hObject, eventdata, handles);
end




% function radiobutton1_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton1,'Value')==0
%     set(handles.radiobutton1,'Value',1);
% else
%     set(handles.radiobutton2,'Value',0);
%     set(handles.radiobutton3,'Value',0);
%     set(handles.lobatto,'Value',0);
% end



function radiobutton2_Callback(hObject, eventdata, handles)

if get(handles.radiobutton2,'Value')==0
    set(handles.radiobutton2,'Value',1);
else
    set(handles.radiobutton1,'Value',0);
    set(handles.radiobutton3,'Value',0);
    set(handles.lobatto,'Value',0);
end


function radiobutton3_Callback(hObject, eventdata, handles)
if get(handles.radiobutton3,'Value')==0
    set(handles.radiobutton3,'Value',1);
else
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton1,'Value',0);
    set(handles.lobatto,'Value',0);
end

% --- Executes on button press in lobatto.
function lobatto_Callback(hObject, eventdata, handles)
if get(handles.lobatto,'Value')==0
    set(handles.lobatto,'Value',1);
else
    set(handles.radiobutton3,'Value',0);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton1,'Value',0);
end

function loeschen_Callback(hObject, eventdata, handles)

      set(handles.edit1,'String','');
      set(handles.checkbox4,'Value',0);
      set(handles.edit2,'String','1');
      set(handles.edit11,'String','');
      set(handles.radiobutton1,'Value',1);
      set(handles.radiobutton2,'Value',0);
      set(handles.radiobutton3,'Value',0);
      set(handles.edit13,'String','');
      set(handles.edit14,'String','');
      set(handles.edit15,'String','');
      set(handles.edit20,'String','');
      set(handles.parameter,'String','');
      set(handles.variablen,'String','');
      set(handles.radiobutton8,'Value',0);
      set(handles.radiobutton12,'Value',0);
      set(handles.edit24,'String','');
      set(handles.c,'String','');
      set([handles.ausfuehren,handles.lobatto,handles.radiobutton12,handles.edit21,handles.pushbutton29,handles.graphischdarstellen,handles.ergebnisse,handles.fehler,handles.pathfollowing,handles.meshadaptation,handles.checkbox4,handles.startprofildatei,handles.startprofillastmat,handles.checkbox2,handles.radiobutton4,handles.einstellungen,handles.popupmenu1,handles.startprofilberechnen,handles.edit2,handles.neuesgitter,handles.startprofilparameter,handles.startprofilwerte,handles.startprofilgitter,handles.edit11,handles.radiobutton1,handles.radiobutton2,handles.radiobutton3,handles.edit13,handles.edit14,handles.edit15,handles.edit20,handles.parameter,handles.variablen,handles.radiobutton8,handles.radiobutton12, handles.edit24],'Enable','on')
      set([handles.startprofilgitter,handles.startprofillastmat,handles.startprofilwerte,handles.neuesgitter,handles.popupmenu1,handles.startprofilberechnen,handles.checkbox2,handles.radiobutton4,],'Enable','on');
    
   set(handles.automatic,'Value',1);
   set(handles.edit13,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   set(handles.edit13,'Enable','off');
  %set(handles.speichern,'Enable','on');
   load options.mat

if isempty(get(handles.edit1,'String'))
 
abstol =1e-10;
reltol = 1e-10;
maxiter= 90000000;
maxfunevals=90000000;
updatejacfactor = 0.5;
lambdamin = 0.001;
switchtoffnfactor = 0.5;
ausgabe =1;
TRM=0;
feinesgitter = 1;
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
   
save options bvpopt abstol reltol maxiter maxfunevals updatejacfactor lambdamin switchtoffnfactor ausgabe TRM abstolgitter reltolgitter K wiederholung feinesgitter n0;
end    
   
 
   if  abstol <= 10^-6

       set(handles.edit13,'String',8)

   elseif abstol > 10^-6 && abstol <= 10^-4

       set(handles.edit13,'String',6);

   elseif abstol > 10^-4 && abstol <= 10^-2

       set(handles.edit13,'String',4);
   elseif abstol > 10^-2

       set(handles.edit13,'String',2);
   end


      
   
     
function ausfuehren_Callback(hObject, eventdata,handles)
meshadaptation=get(handles.meshadaptation,'Value');
%msgbox('Wait a moment, ...','Info','modal');
if meshadaptation
    auswahl=3;
else
 
    auswahl=1;
  
end

if get(handles.radiobutton8,'Value')==1 
    antwort_EVP='true';
  

end




if get(handles.radiobutton12,'Value')==1
    antwort_Infinite='true';
end 


if get(handles.automatic,'Value')==1

    antwort_auto='true';
    
    
else
     antwort_auto='false';
    
end

try
    [coeff,x1,valx1,x1tau,valx1tau,polynomials,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,error,lambda,eigenfunction] = berechnen(auswahl,hObject,eventdata,handles);
catch e
    fprintf(1,'%s\n',e.message);
    msgbox(e.message,'Error in solver','error');
    return;
end
    
main_gui_handle=bvpsuite;
mainGUIdata  = guidata(main_gui_handle);



pfad=get(mainGUIdata.pfad,'String');

aktuellesverzeichnis=cd;
arbeitsverzeichnis=pfad;
cd(arbeitsverzeichnis);


filename=get(mainGUIdata.edit1,'String');
bvpfile=strrep(filename,'.m','');


parameter=str2num(get(handles.parameter,'String'));

if parameter>0
    p=coeff(length(coeff)-parameter+1:length(coeff));
    p=p';
    if length(strfind(char(version),'R14'))==0
        
        if feval(bvpfile,'Infinite')


            if feval(bvpfile,'EVP')

                cd(aktuellesverzeichnis);
                save last coeff x1 valx1 x1tau valx1tau polynomials p  sol_infinite tau_infinite x1_infinite solx1_infinite lambda eigenfunction;
            else

                cd(aktuellesverzeichnis)
                save last coeff x1 valx1 x1tau valx1tau polynomials p  sol_infinite tau_infinite solx1_infinite x1_infinite;
            end

        else

            if feval(bvpfile,'EVP')
                cd(aktuellesverzeichnis)
                save last coeff x1 valx1 x1tau valx1tau polynomials p lambda eigenfunction;
            
            else
                cd(aktuellesverzeichnis)
                save last coeff x1 valx1 x1tau valx1tau polynomials p;

            end

        end
        
    else
        
        if feval(bvpfile,'Infinite')


            if feval(bvpfile,'EVP')
                    
                cd(aktuellesverzeichnis)
                save last -v6 coeff x1 valx1 x1tau valx1tau polynomials p  sol_infinite tau_infinite solx1_infinite x1_infinite  lambda eigenfunction;
            else

                cd(aktuellesverzeichnis)
                save last -v6 coeff x1 valx1 x1tau valx1tau polynomials p  sol_infinite tau_infinite solx1_infinite x1_infinite;
            end

        else

            if feval(bvpfile,'EVP')
                
                cd(aktuellesverzeichnis)
                save last -v6  coeff x1 valx1 x1tau valx1tau polynomials p lambda eigenfunction;
            
            else
                cd(aktuellesverzeichnis)
                save last -v6  coeff x1 valx1 x1tau valx1tau polynomials p;

            end

        end
    
    end
    if (auswahl==1 || auswahl==3)
        if error==1
      
        else    
               helpdlg('The calculation was successful. The values are saved in last.mat. They are now available in the workspace. For further information see "Workspace" in the Matlab help.','SUCCESS');
        end
    end
else
    if length(strfind(char(version),'R14'))==0
        
        if feval(bvpfile,'Infinite')


            if feval(bvpfile,'EVP')
                cd(aktuellesverzeichnis)
                 save last coeff x1 valx1 x1tau valx1tau polynomials sol_infinite tau_infinite solx1_infinite x1_infinite lambda eigenfunction ;
            else
                
                cd(aktuellesverzeichnis)
                save last coeff x1 valx1 x1tau valx1tau polynomials sol_infinite tau_infinite  solx1_infinite x1_infinite;
            end

        else

            if feval(bvpfile,'EVP')
                cd(aktuellesverzeichnis)
                save last coeff x1 valx1 x1tau valx1tau polynomials  lambda eigenfunction;
            
            else
                cd(aktuellesverzeichnis)
                save last coeff x1 valx1 x1tau valx1tau polynomials ;

            end

        end
        
    else
        
        if feval(bvpfile,'Infinite')


            if feval(bvpfile,'EVP')

                cd(aktuellesverzeichnis)
                save last -v6 coeff x1 valx1 x1tau valx1tau polynomials sol_infinite tau_infinite solx1_infinite x1_infinite  lambda eigenfunction;
            else

                cd(aktuellesverzeichnis)
                save last -v6 coeff x1 valx1 x1tau valx1tau polynomials sol_infinite tau_infinite solx1_infinite x1_infinite ;
            end

        else

            if feval(bvpfile,'EVP')
                
                cd(aktuellesverzeichnis)
                save last -v6  coeff x1 valx1 x1tau valx1tau polynomials lambda eigenfunction;
            
            else
                
                cd(aktuellesverzeichnis)
                save last -v6  coeff x1 valx1 x1tau valx1tau polynomials;

            end

        end
    
    end
    if (auswahl==1 || auswahl==3)
        
        if error==1
      
        else    
               helpdlg('The calculation was successful. The values are saved in last.mat. They are now available in the workspace. For further information see "Workspace" in the Matlab help.','SUCCESS');
        end
        
    end
end

function konvergenzord_Callback(hObject, eventdata, handles)
auswahl=2;


berechnen(auswahl,hObject,eventdata,handles);

  



function [coeff,x1,valx1,x1tau,valx1tau,polynomials,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,error,lambda,eigenfunction]=berechnen(auswahl,hObject,eventdata,handles);
%Die Routine wird sowohl f�r das Berechnen der L�sung �ber einem speziellen
%Gitter als auch f�r die Initialisierung der Konvergenzordnung oder des Pathfollowing verwendet.
eigenfunction=[];
error=0;
lambda=[];
sol_infinite=[];
tau_infinite=[];
if get(handles.zeichnen,'Value')==1
    zeichnen=1;
else
    zeichnen=0;
end
bvpfile=get(handles.edit1,'String');
bvpfilem=bvpfile;
bvpfile=strrep(bvpfile,'.m','');
pfad=get(handles.pfad,'String');
% Es gibt ein Problem mit den Function-Handles: Wird innerhalb von
% bvpsuite eine Datei �berschrieben, die bereits existiert, so bleiben
% alle Function-Handles, die schon einmal von der urspr�nglichen Datei
% ausgef�hrt wurden, erhalten - werden also nicht durch die neue Datei
% ersetzt. Das h�tte zur Folge, da� man vor dem Dr�cken von "Berechnen",
% das ganze Programm einmal schlie�en und wieder �ffnen m��te, damit die
% neuen Daten, die mit "speichern" gespeichert wurden, auch bei der
% Berechnung aktiv werden. Um dieses Problem zu umgehen, wird eine Datei
% mit einem zuf�lligen Namen (bestehend aus Datum und einer 5stelligen
% Zufallszahl - eine Kollission der Namen ist ausgeschlossen) erstellt, einmal mit feval ausgewertet, um das
% Function-Handle zu speichern, und dann sofort wieder gel�scht. Das
% Function-Handle bleibt aber erhalten und kann weiterverwendet werden.

%Aufgrund der neuen Speicherstruktur ist es nicht mehr notwendig, einen
%Zufallsnamen zu generieren.
%zufallsname=strcat('temp_bvpsuite_file_',strcat(strcat(strrep(strrep(num2str(clock),' ',''),'.',''),strrep(num2str(randsrc(1,5,[0 1 2 3 4 5 6 7 8 9])),' ','')),'.m'));
%zufallsname='temp_bvpsuite_file.m';
datei=fopen(strcat(pfad,bvpfilem),'r');
if datei==-1
    err('bvps_errdlg17');
    return;
end
inhaltdatei=fread(datei);
fclose(datei);
inhaltdateibak=inhaltdatei;
inhaltdatei=char(inhaltdatei');
if (length(strfind(inhaltdatei,'Kontrollnummer43753976430976655143'))>0)
    err('bvps_errdlg18');
    return;
end
gueltig=strfind(inhaltdatei,'Kontrollnummer43753976430976655144');
if (length(gueltig)>0) 
      err('bvps_errdlg18');
      return;
end

if (length(gueltig)==0 & length(strfind(inhaltdatei,'Kontrollnummer43753976430976655145'))==0  && length(strfind(inhaltdatei,'Kontrollnummer43753976430976655143'))==0)
      errordlg('The chosen file is not compatible with this version of bvpsuite!','Error');
      return;
end
checkm3=strfind(bvpfilem,'bvpsuite.m');
if(length(checkm3)~=0)
    errordlg('bvpsuite is not a bvpfile!','Error');
    return;
end
%delete('temp_bvpsuite_file_*.m');
%datei=fopen(zufallsname,'w');
%fprintf(datei,'%s',char(inhaltdateibak));
%fclose(datei);
%copyfile mu�te wegen Komplikationen mit Unix entfernt werden
%copyfile(bvpfilem,zufallsname);
%zufallsname2=strrep(zufallsname,'.m','');
%Dummyauswertung
addpath(pfad);
feval(bvpfile,'x1');
%L�schen (Function-Handle bleibt wegen feval erhalten)
%delete(zufallsname);
load options;
if get(handles.radiobutton4,'Value')==1 %radiobutton4 ist "Gespeichertes Profil"
    if length(strfind(inhaltdatei,'startprofilgitter'))~=0 %Ein Startprofil ist in der Datei gespeichert
        help=inhaltdatei(strfind(inhaltdatei,'%startprofilgitter'):strfind(inhaltdatei,'#startprofilgitter')-1);
        help=strrep(help,'%startprofilgitter','');
        help=help(3:length(help));
        startprofilgitter=str2num(help);
        help=inhaltdatei(strfind(inhaltdatei,'%startprofilwerte'):strfind(inhaltdatei,'#startprofilwerte')-1);
        help=strrep(help,'%startprofilwerte','');
        help=help(3:length(help));
        startprofilwerte=str2num(help);
        if length(strfind(inhaltdatei,'startprofilparameter'))~=0
             help=inhaltdatei(strfind(inhaltdatei,'%startprofilparameter'):strfind(inhaltdatei,'%#startprofilparameter')-1);
             help=strrep(help,'%startprofilparameter','');
             help=help(3:length(help));
             startprofilparameter=str2num(help);
        else
             startprofilparameter=[];
        end
        
        if length(strfind(inhaltdatei,'lambda'))~=0
             help=inhaltdatei(strfind(inhaltdatei,'%lambda'):strfind(inhaltdatei,'%#lambda')-1);
             help=strrep(help,'%lambda','');
             help=help(3:length(help));
             startprofillambda=str2num(help);
        else
             startprofillambda=[];
        end    
       
        load initialmesh.mat;
        
        %if min(size(startprofilwerte)==size(werte))==0 dimension of values matrices do not agree
        %if  length(startprofilgitter)~=length(stellen) dimension of mesh
        %vectors do not agree
        %if length(feval(bvpfile,'x1'))~=length(x1) size of calculation mesh
        %is different
        %if feval(bvpfile,'standard')~=standard  if not standard
        %coll.points 
        %if length(feval(bvpfile,'rho'))~=length(rho) if different
        %coll. points
        %if length(num2str(startprofilparameter))~=length(num2str(p)) if
        %number of parameters is different 
        
       
        %if the sizes are not the same 
        
        
        
        
        
        
        
%         if min(size(startprofilwerte)==size(werte))==0 || length(startprofilgitter)~=length(stellen) || length(feval(bvpfile,'x1'))~=length(x1) || feval(bvpfile,'standard')~=standard || length(feval(bvpfile,'rho'))~=length(rho) || length(num2str(startprofilparameter))~=length(num2str(p)) || length(num2str(startprofillambda))~=length(num2str(lambda))
%             msgbox('Recalculate initial mesh ...','Info','modal');
%             
%             
%             if feval(bvpfile,'Infinite')
%             
%              help=inhaltdatei(strfind(inhaltdatei,'%x1'):strfind(inhaltdatei,'%#x1')-1);
%              help=strrep(help,'%x1','');
%              help=help(3:length(help));
%              x1=str2num(help);
%               [x1,start]=initialmesh(bvpfile,startprofilgitter,startprofilwerte,x1,startprofilparameter,startprofillambda);
%                 
%             else    
%             [x1,start]=initialmesh(bvpfile,startprofilgitter,startprofilwerte,feval(bvpfile,'x1'),startprofilparameter,startprofillambda);
%                       
%             
%             end 
%             
%             stellen=startprofilgitter;
%             werte=startprofilwerte;
%             standard=feval(bvpfile,'standard');
%             rho=feval(bvpfile,'rho');
%             p=startprofilparameter;
%             lambda=startprofillambda;
%             if length(strfind(char(version),'R14'))==0
%                 save initialmesh x1 start stellen werte standard rho p lambda;
%             else
%                 save initialmesh -v6 x1 start stellen werte standard rho p lambda;
%            end
            
            
            
      %  elseif min(min(startprofilwerte==werte))==0 || min(startprofilgitter==stellen)==0 || min(feval(bvpfile,'x1')==x1)==0 || min(feval(bvpfile,'rho')==rho)==0 || min(strcat('[',num2str(startprofilparameter))==strcat('[',num2str(p)))==0 || strcmp(num2str(startprofillambda),num2str(lambda))==0
%         elseif strcmp(num2str(startprofilwerte),num2str(werte))==0 || strcmp(num2str(startprofilgitter),num2str(stellen))==0  || min(feval(bvpfile,'x1')==x1)==0 || min(feval(bvpfile,'rho')==rho)==0 || min(strcat('[',num2str(startprofilparameter))==strcat('[',num2str(p)))==0 || strcmp(num2str(startprofillambda),num2str(lambda))==0
         %if sizes are the same but different values
         %min(min(startprofilwerte==werte))==0 %values are not the same
         %min(startprofilgitter==stellen)==0   %mesh is not the same
         %min(feval(bvpfile,'x1')==x1)==0 %x1 is different
         %min(feval(bvpfile,'rho')==rho)==0 
         %min(strcat('[',num2str(startprofilparameter))==strcat('[',num2str(p)))==0   
            
            
            
            %msgbox('Recalculate initial mesh ...','Info','modal');
            
            if feval(bvpfile,'Infinite')
            
             help=inhaltdatei(strfind(inhaltdatei,'%x1'):strfind(inhaltdatei,'%#x1')-1);
             help=strrep(help,'%x1','');
             help=help(3:length(help));
             x1=str2num(help);
                
            [x1,start]=initialmesh(bvpfile,startprofilgitter,startprofilwerte,x1,startprofilparameter,startprofillambda,0);
                
            else 
            [x1,start]=initialmesh(bvpfile,startprofilgitter,startprofilwerte,feval(bvpfile,'x1'),startprofilparameter,startprofillambda,0);

            end 
 
            stellen=startprofilgitter;
            werte=startprofilwerte;
            p=startprofilparameter;
            standard=feval(bvpfile,'standard');
            rho=feval(bvpfile,'rho');
            %Dasha EDIT- just take first Lambda of array (dont loop
            %through)
            lambda=startprofillambda;
            if length(strfind(char(version),'R14'))==0
                save initialmesh x1 start stellen werte standard rho p lambda;
            else
                save initialmesh -v6 x1 start stellen werte standard rho p lambda;
            end
        
  %      end
        
      
      
        %if none of the two cases-> then take initialmesh.mat 
        switch auswahl
            case 1
                [coeff,x1,valx1,x1tau,valx1tau,polynomials,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,lambda,eigenfunction]=run(bvpfile,zeichnen,x1,start,bvpopt,ausgabe);
              
               
                
                
            case 2
                if length(strfind(char(version),'R14'))==0
                    save help1 bvpfile zeichnen x1 start bvpopt ausgabe;
                else
                    save help1 -v6 bvpfile zeichnen x1 start bvpopt ausgabe;
                end
                convergencegui;
                coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;
            case 3                
                [coeff,x1,valx1,x1tau,valx1tau,polynomials,valerror,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,lambda,eigenfunction]=meshadaptation(abstolgitter,reltolgitter,K,bvpfile,zeichnen,x1,start,bvpopt,ausgabe,wiederholung,feinesgitter);
            case 4
                if length(strfind(char(version),'R14'))==0
                    save help1 bvpfile zeichnen bvpopt ausgabe;
                else
                    save help1 -v6 bvpfile zeichnen bvpopt ausgabe;
                end
                pathfollowinggui;
                coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;
            
        end
    else %wenn kein Startprofil in der Datei gespeichert ist
        switch auswahl
            case 1
                [coeff,x1,valx1,x1tau,valx1tau,polynomials,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,lambda,eigenfunction]=run(bvpfile,zeichnen,[],[],bvpopt,ausgabe);
                 
              
          
            case 2
                x1=[];
                start=[];
                if length(strfind(char(version),'R14'))==0
                    save help1 bvpfile zeichnen x1 start bvpopt ausgabe;
                else
                    save help1 -v6 bvpfile zeichnen x1 start bvpopt ausgabe;
                end
                convergencegui;
                coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;
            case 3
                 [coeff,x1,valx1,x1tau,valx1tau,polynomials,valerror,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,lambda,eigenfunction]=meshadaptation(abstolgitter,reltolgitter,K,bvpfile,zeichnen,[],[],bvpopt,ausgabe,wiederholung,feinesgitter);
            case 4
                if length(strfind(char(version),'R14'))==0
                    save help1 bvpfile zeichnen bvpopt ausgabe;
                else
                    save help1 -v6 bvpfile zeichnen bvpopt ausgabe;
                end
                pathfollowinggui;
                coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;
            
        end
    end
end

%Profile in the Fields!!
if get(handles.checkbox2,'Value')==1 
        
     
    if feval(bvpfile,'EVP')
      
         if length(get(handles.lambda,'String'))==0
            err('bvps_errdlg20');
            return;
        end
    
    else
        if length(get(handles.startprofilgitter,'String'))==0
            err('bvps_errdlg21');
            return;
        end
        if length(get(handles.startprofilwerte,'String'))==0
            err('bvps_errdlg22');
            return;
        end

        if length(get(handles.neuesgitter,'String'))==0
            err('bvps_errdlg23');
            return;
        end

    end     
        

        load initialmesh.mat;
        neuestellen=str2num(get(handles.startprofilgitter,'String'));
        neuewerte=str2num(get(handles.startprofilwerte,'String'));
        
        if feval(bvpfile,'EVP') && length(get(handles.startprofilgitter,'String'))==0
            neuesx1=feval(bvpfile,'x1');
        else
         neuesx1=str2num(get(handles.neuesgitter,'String'));    
        end
        neuesp=str2num(get(handles.startprofilparameter,'String'));
        neueslambda=str2num(get(handles.lambda,'String'));


%         if min(size(neuewerte)==size(werte))==0 || length(neuestellen)~=length(stellen) || length(neuesx1)~=length(x1) || feval(bvpfile,'standard')~=standard || length(feval(bvpfile,'rho'))~=length(rho) || length(num2str(p))~=length(num2str(neuesp)) || length(num2str(neueslambda))~=length(num2str(lambda))
%             msgbox('Recalculation of the initial profile ...','Info','modal');
% 
% 
% 
% 
%              [x1,start]=initialmesh(bvpfile,neuestellen,neuewerte,neuesx1,neuesp,neueslambda);
% 
% 
%             stellen=neuestellen;
%             werte=neuewerte;
%             standard=feval(bvpfile,'standard');
%             rho=feval(bvpfile,'rho');
%             p=neuesp;
%             lambda=neueslambda;
%             if length(strfind(char(version),'R14'))==0
%                 save initialmesh x1 start stellen werte standard rho p lambda;
%             else
%                 save initialmesh -v6 x1 start stellen werte standard rho p lambda ;
%             end
%         %elseif min(min(neuewerte==werte))==0 || min(neuestellen==stellen)==0 || min(neuesx1==x1)==0 || min(feval(bvpfile,'rho')==rho)==0 || min(strcat('[',num2str(neuesp))==strcat('[',num2str(p)))==0 || strcmp(num2str(neueslambda),num2str(lambda))==0
%          elseif strcmp(num2str(neuewerte),num2str(werte))==0 || strcmp(num2str(neuestellen),num2str(stellen))==0 ||strcmp(num2str(neuesx1),num2str(x1))==0 || min(feval(bvpfile,'rho')==rho)==0 || min(strcat('[',num2str(neuesp))==strcat('[',num2str(p)))==0 || strcmp(num2str(neueslambda),num2str(lambda))==0   
            %msgbox('Recalculation of the initial profile ...','Info','modal');
% Dasha EDIT FOR LOOP WHEN Lambda is an array
ErrorArray = zeros(1,size(neueslambda,2));
LambdaArray = zeros(1,size(neueslambda,2));
if isempty(LambdaArray)
    ArraySize = 1;
else
    ArraySize = size(neueslambda,2);
end
    for lambdaCount = 1:ArraySize
            
            if isempty(LambdaArray)
                lambdaVal = neueslambda;
            else
                lambdaVal = neueslambda(lambdaCount);
            end
            [x1,start]=initialmesh(bvpfile,neuestellen,neuewerte,neuesx1,neuesp,lambdaVal,0);

            stellen=neuestellen;
            werte=neuewerte;
            standard=feval(bvpfile,'standard');
            rho=feval(bvpfile,'rho');
            p=neuesp;
            if isempty(LambdaArray)
                lambda=neueslambda;
            else
                lambda=neueslambda(lambdaCount);
            end
            fprintf('Eigenvalue initial = %8e\n',lambda);
            if length(strfind(char(version),'R14'))==0
                save initialmesh x1 start stellen werte standard rho p lambda;
            else
                save initialmesh -v6 x1 start stellen werte standard rho p lambda;
            end

%        end


 
        valerror = 1;
        switch auswahl
            case 1


                [coeff,x1,valx1,x1tau,valx1tau,polynomials,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,Lambda,eigenfunction]=run(bvpfile,zeichnen,x1,start,bvpopt,ausgabe);


            case 2
                if length(strfind(char(version),'R14'))==0
                    save help1 bvpfile zeichnen x1 start bvpopt ausgabe;
                else
                    save help1 -v6 bvpfile zeichnen x1 start bvpopt ausgabe;
                end
                %convergencegui;
                coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;
            case 3
                [coeff,x1,valx1,x1tau,valx1tau,polynomials,valerror,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,Lambda,eigenfunction]=meshadaptation(abstolgitter,reltolgitter,K,bvpfile,zeichnen,x1,start,bvpopt,ausgabe,wiederholung,feinesgitter);            
            case 4
                if length(strfind(char(version),'R14'))==0
                    save help1 bvpfile zeichnen bvpopt ausgabe;
                else
                    save help1 -v6 bvpfile zeichnen bvpopt ausgabe;
                end
                pathfollowinggui;
                coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;

        end
        
%         if length(polynomials) > 0
% 
%             %Dasha Edit (save these files separately from last.mat)
%             %EVPfileName = strcat('lambda',num2str(lambda),'to',num2str(valx1tau(2,1)),'.mat');
             EVPfileName = strcat('lambda',num2str(lambda),'.mat');
             save(EVPfileName,'valx1tau','x1tau');
%             ErrorArray(lambdaCount) = polynomials(1);
%         else
%             ErrorArray(lambdaCount) = -1;
%         end
        %LambdaArray(lambdaCount) = lambda;
        %printf('Eigenvalue initial = %8e\n',lambda);
    end
    %save ErrorLambdas ErrorArray LambdaArray
end



if get(handles.startprofillastmat,'Value')==1 %As Profile 
 
  
    pfad=char(get(handles.pfad,'String'));
    coeff=[];
    start=[];
    p=[];
    startprofildatei=get(handles.startprofildatei,'String');
 
    if length(startprofildatei)~=4 || ~min(startprofildatei=='last')
  
        load(strcat(pfad,startprofildatei))
    else
        load(startprofildatei)
    end
    
  
    if length(start)>0
        coeff=start;
    end
    if length(get(handles.neuesgitter,'String'))==0
       
         error=1;       
         err('bvps_errdlg23');
        return;
    end
    neuesx1=str2num(get(handles.neuesgitter,'String'));
   
  
    if length(p)>0
        
        if feval(bvpfile,'Infinite') 
           
            if feval(bvpfile,'EVP')
            
            else
            [x1,start]=initialmesh(bvpfile,tau_infinite,sol_infinite,neuesx1,p,lambda,0);
        
            end 
        else
           
            
        if feval(bvpfile,'EVP')   
            
            [x1,start]=initialmesh(bvpfile,x1tau,valx1tau,neuesx1,p,lambda,0);
        else    
            [x1,start]=initialmesh(bvpfile,x1tau,valx1tau,neuesx1,p,0);
        end 
        
        
        end
    else
        
        if feval(bvpfile,'Infinite')  
      
            if feval(bvpfile,'EVP')
            
                 [x1,start]=initialmesh(bvpfile,tau_infinite,sol_infinite,neuesx1,[],lambda,0);
            else
                [x1,start]=initialmesh(bvpfile,tau_infinite,sol_infinite,neuesx1,[],[],0);
            end
        else

            if feval(bvpfile,'EVP')
                              
                            
                [x1,start]=initialmesh(bvpfile,x1tau,valx1tau,neuesx1,[],lambda,0);
            else    
                [x1,start]=initialmesh(bvpfile,x1tau,valx1tau,neuesx1,[],[],0);
            end 
            
        end 
    end
    stellen=0;
    werte=0;
    
    standard=feval(bvpfile,'standard');
    rho=feval(bvpfile,'rho');
    p_n=feval(bvpfile,'parameter');
%    p=coeff(end-p_n+1:end);
    
    if length(strfind(char(version),'R14'))==0
        save initialmesh x1 start stellen werte standard rho p lambda;
    else
        save initialmesh -v6 x1 start stellen werte standard rho p lambda;
    end
    switch auswahl
        case 1
            [coeff,x1,valx1,x1tau,valx1tau,polynomials,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,lambda,eigenfunction]=run(bvpfile,zeichnen,x1,start,bvpopt,ausgabe);
        
 
        
        
        case 2
            if length(strfind(char(version),'R14'))==0
                save help1 bvpfile zeichnen x1 start bvpopt ausgabe;
            else
                save help1 -v6 bvpfile zeichnen x1 start bvpopt ausgabe;
            end
            convergencegui;
            coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;
        case 3
            [coeff,x1,valx1,x1tau,valx1tau,polynomials,valerror,sol_infinite,tau_infinite,solx1_infinite,x1_infinite,lambda,eigenfunction]=meshadaptation(abstolgitter,reltolgitter,K,bvpfile,zeichnen,x1,start,bvpopt,ausgabe,wiederholung,feinesgitter);            
        case 4
            if length(strfind(char(version),'R14'))==0
                save help1 bvpfile zeichnen bvpopt ausgabe;
            else
                save help1 -v6 bvpfile zeichnen bvpopt ausgabe;
            end
            pathfollowinggui;
            coeff=0;x1=0;valx1=0;x1tau=0;valx1tau=0;polynomials=0;

    end


end








function zeichnen_Callback(hObject, eventdata, handles)


% --- Executes on button press in ergebnisse.
function ergebnisse_Callback(hObject, eventdata, handles)
    datei=fopen('last.mat');
    if datei==-1
        err('bvps_errdlg24');
        return;
    end
    fclose(datei);
load last.mat;
parameter=str2num(get(handles.parameter,'String'));
if parameter>0
    p=coeff(length(coeff)-parameter+1:length(coeff));
    p=p';
    assignin('base','parameter',p);
end


if get(handles.radiobutton12,'Value') & get(handles.radiobutton8,'Value')==0
    assignin('base','sol_infinite',sol_infinite);
    assignin('base','tau_infinite',tau_infinite);
assignin('base','solx1_infinite',solx1_infinite);
    assignin('base','x1_infinite',x1_infinite);

elseif get(handles.radiobutton12,'Value') & get(handles.radiobutton8,'Value')==1
    assignin('base','tau_infinite',tau_infinite);
  assignin('base','lambda',lambda);
 assignin('base','eigenfunction',eigenfunction);  
elseif get(handles.radiobutton12,'Value') ==0 & get(handles.radiobutton8,'Value')==1    
 assignin('base','lambda',lambda);
 assignin('base','eigenfunction',eigenfunction);

else
    
    
assignin('base','coeff',coeff);
assignin('base','x1',x1);
assignin('base','valx1',valx1);
assignin('base','x1tau',x1tau);
assignin('base','valx1tau',valx1tau);
%assignin('base','polynomials',polynomials);

    
end
 

load options.mat;

%Das "parameter" in der if Abfrage ist jenes der Datei (Skalar), das andere (open) jenes
%des Workspace (i.d.R. ein Vektor)
if parameter>0
    open parameter;
end



if get(handles.radiobutton12,'Value') & get(handles.radiobutton8,'Value')==0
     open solx1_infinite;
    open x1_infinite;
    open sol_infinite;
    open tau_infinite;
elseif get(handles.radiobutton12,'Value') & get(handles.radiobutton8,'Value')==1
    open tau_infinite;
    open lambda;
    open eigenfunction;    
elseif get(handles.radiobutton12,'Value') ==0 & get(handles.radiobutton8,'Value')==1    
    open valx1tau;
    open lambda;
    open eigenfunction;    
else 

 open x1;
%open valx1;
open x1tau;
open valx1tau;
%open polynomials; 
    
end 

function startprofilgitter_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function startprofilgitter_Callback(hObject, eventdata, handles)


function startprofilwerte_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function startprofilwerte_Callback(hObject, eventdata, handles)


function neuesgitter_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function neuesgitter_Callback(hObject, eventdata, handles)


function checkbox2_Callback(hObject, eventdata, handles)
set(handles.radiobutton4,'value',0);
set(handles.checkbox2,'value',1);
set(handles.startprofillastmat,'value',0);

function radiobutton4_Callback(hObject, eventdata, handles)
set(handles.checkbox2,'value',0);
set(handles.radiobutton4,'value',1);
set(handles.startprofillastmat,'value',0);

function startprofillastmat_Callback(hObject, eventdata, handles)
set(handles.checkbox2,'value',0);
set(handles.radiobutton4,'value',0);
set(handles.startprofillastmat,'value',1);

function startprofilberechnen_Callback(hObject, eventdata, handles)
    bvpfile=get(handles.edit1,'String');
    pfad=get(handles.pfad,'String');
    bvpfilem=bvpfile;
    bvpfile=strrep(bvpfile,'.m','');
    datei=fopen(strcat(pfad,bvpfilem),'r');
    if datei==-1
        err('bvps_errdlg17');
        return;
    end
    inhaltdateibak=fread(datei);
    inhaltdatei=char(inhaltdateibak');
    fclose(datei);    
%zufallsname=strcat('temp_bvpsuite_file_',strcat(strcat(strrep(strrep(num2str(clock),' ',''),'.',''),strrep(num2str(randsrc(1,5,[0 1 2 3 4 5 6 7 8 9])),' ','')),'.m'));
%zufallsname='temp_bvpsuite_file.m';
%delete('temp_bvpsuite_file_*.m');
%datei=fopen(zufallsname,'w');
%fprintf(datei,'%s',char(inhaltdateibak));
%fclose(datei);
%copyfile(bvpfilem,zufallsname);
%zufallsname2=strrep(zufallsname,'.m','');
%Dummyauswertung
feval(bvpfile,'x1');
%L�schen (Function-Handle bleibt wegen feval erhalten)
%delete(zufallsname);

%handles.popupmenu1 defines if 

%1.show initial values in fields above
%2.show saved initial values
%3. use last solution


if get(handles.popupmenu1,'Value')==1 %1.show initial values in fields above
    datei=fopen('initialmesh.mat');
    if datei==-1
        errordlg('The file "initialmesh.mat" does not exist, calculate an initial profile first!');
        return;
    end
    fclose(datei);
    if length(get(handles.startprofilgitter,'String'))==0
        err('bvps_errdlg21');
        return;
    end
    if length(get(handles.startprofilwerte,'String'))==0
        err('bvps_errdlg22');
        return;
    end
    if length(get(handles.neuesgitter,'String'))==0
        err('bvps_errdlg23');
        return;
    end
    load initialmesh.mat;
    

   
    
    %reads initial profile from GUI
    neuestellen=str2num(get(handles.startprofilgitter,'String'));
    neuewerte=str2num(get(handles.startprofilwerte,'String'));
    neuesx1=str2num(get(handles.neuesgitter,'String')); %alternative mesh in GUI
    neuesp=str2num(get(handles.startprofilparameter,'String'));
    neueslambda=str2num(get(handles.lambda,'String'));
    %Dasha Edit
    if isempty(neueslambda)
        ArraySize = 1;
    else
        ArraySize = size(neueslambda,2);
    end
    for lambdaCount = 1:ArraySize
        if isempty(neueslambda)
            lambda=neueslambda;
        else
            lambda=neueslambda(lambdaCount);
        end
        
        [x1,start]=initialmesh(bvpfile,neuestellen,neuewerte,neuesx1,neuesp,lambda,0);



        stellen=neuestellen;
        werte=neuewerte;
        standard=feval(bvpfile,'standard');
        rho=feval(bvpfile,'rho');
        p=neuesp;
        
        if length(strfind(char(version),'R14'))==0
            save initialmesh x1 start stellen werte standard rho p lambda;
        else
            save initialmesh -v6 x1 start stellen werte standard rho p lambda;
        end
        %    end
        % figure


        %polynomials=equations('polynom',start,bvpfile,x1);
        valx1tau=equations('valx1tau',start,bvpfile,x1);
        x1tau=equations('x1tau',start,bvpfile,x1);
        %plotpoly(polynomials,x1,feval(bvpfile,'n'));
        n=feval(bvpfile,'n');
  %  end
    
    
    
    
    
   
    ep=str2num(get(handles.edit24,'String')) ;
   
    figure
    
    if feval(bvpfile,'Infinite')&& feval(bvpfile,'EVP')==0


        if ep ~=0

            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end

            for i=1:n
                
                subplot(n,1,i);
                plot(tau_infinite,sol_infinite(i,:),'k');
                xlim([stellen(1),stellen(end)])
                
                 xlabel('t');
                 ylabel(['z_',num2str(i)]);
            end

            title(subplot(n,1,1),'Initial values');
        else

            for i=1:n/2


                
                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),0);

            end

            for i=1:n/2

                subplot(n/2,1,i);


                plot(tau_infinite,sol_infinite(i,:),'k');
                xlim([stellen(1),stellen(end)])
                
                 xlabel('t');
                 ylabel(['z_',num2str(i)]);

            end

            
            title(subplot(n/2,1,1),'Initial values');
            
        end

    elseif feval(bvpfile,'Infinite') && feval(bvpfile,'EVP')



        if ep == 0


            for i=1:n/2



                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);



            end


            for i=1:n/2-2
                subplot(n/2-2,1,i);
                plot(tau_infinite,sol_infinite(i,:),'k');
                xlim([stellen(1),stellen(end)])
               
                a=sol_infinite(n/2-1,1);
                xlabel('t');
                ylabel(['z_',num2str(i)]);


            end

            title(subplot(n/2-2,1,1),'Initial values');

        else


            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end


            for i=1:n-2
                subplot(n-2,1,i);



                plot(tau_infinite,sol_infinite(i,:),'k');
                xlim([stellen(1),stellen(end)])
                a=sol_infinite(n-1,1);
                xlabel('t');
                ylabel(['z_',num2str(i)]);

            end

            title(subplot(n-2,1,1),'Initial values');
        end





    elseif feval(bvpfile,'Infinite') ==0 && feval(bvpfile,'EVP')



        for i=1:n-2
            subplot(n-2,1,i);
            plot(x1tau,valx1tau(i,:),'color','black');



            a=valx1tau(n-1,1);
            axis tight
            xlabel('t');
            ylabel(['z_',num2str(i)]);
        end

        title(subplot(n-2,1,1),'Initial values');


    elseif feval(bvpfile,'Infinite') ==0 && feval(bvpfile,'EVP')== 0
        
        for i=1:n
            subplot(n,1,i);
            plot(x1tau,valx1tau(i,:),'color','black');
            axis tight
            xlabel('t');
            ylabel(['z_',num2str(i)]);
        end

title(subplot(n,1,1),'Initial values');
    end
    end


end

if get(handles.popupmenu1,'Value')==2 %2. show saved initial values

    if length(strfind(inhaltdatei,'startprofilgitter'))~=0
        help=inhaltdatei(strfind(inhaltdatei,'%startprofilgitter'):strfind(inhaltdatei,'#startprofilgitter')-1);
        help=strrep(help,'%startprofilgitter','');
        help=help(3:length(help));
        startprofilgitter=str2num(help);
        set(handles.startprofilgitter,'String',help);
        help=inhaltdatei(strfind(inhaltdatei,'%startprofilwerte'):strfind(inhaltdatei,'#startprofilwerte')-1);
        help=strrep(help,'%startprofilwerte','');
        help=help(3:length(help));
        startprofilwerte=str2num(help);
        set(handles.startprofilwerte,'String',help);

       
        if length(strfind(inhaltdatei,'startprofilparameter'))~=0
            help=inhaltdatei(strfind(inhaltdatei,'%startprofilparameter'):strfind(inhaltdatei,'%#startprofilparameter')-1);
            help=strrep(help,'%startprofilparameter','');
            help=help(3:length(help));
            startprofilparameter=str2num(help);
            set(handles.startprofilparameter,'String',help);
        else
            set(handles.startprofilparameter,'String','');
            startprofilparameter=[];
        end


        if length(strfind(inhaltdatei,'lambda'))~=0
            help=inhaltdatei(strfind(inhaltdatei,'%lambda'):strfind(inhaltdatei,'%#lambda')-1);
            help=strrep(help,'%lambda','');
            help=help(3:length(help));
            startprofillambda=str2num(help);
            set(handles.lambda,'String',help);
        else
            set(handles.lambda,'String','');
            startprofillambda=[];
        end


        load initialmesh.mat;

        if feval(bvpfile,'Infinite')

            help=inhaltdatei(strfind(inhaltdatei,'%x1'):strfind(inhaltdatei,'%#x1')-1);
            help=strrep(help,'%x1','');
            help=help(3:length(help));
            x1=str2num(help);
            [x1,start]=initialmesh(bvpfile,startprofilgitter,startprofilwerte,x1,startprofilparameter,startprofillambda,0);

        else
            [x1,start]=initialmesh(bvpfile,startprofilgitter,startprofilwerte,feval(bvpfile,'x1'),startprofilparameter,startprofillambda,0);


        end


        stellen=startprofilgitter;
        werte=startprofilwerte;
        p=startprofilparameter;
        standard=feval(bvpfile,'standard');
        rho=feval(bvpfile,'rho');
        lambda=startprofillambda;
        if length(strfind(char(version),'R14'))==0
            save initialmesh x1 start stellen werte standard rho p lambda;
        else
            save initialmesh -v6 x1 start stellen werte standard rho p lambda;
        end
        % end

    else
        err('bvps_errdlg26');
        return;
    end

    set(handles.checkbox2,'Value',0);
    set(handles.radiobutton4,'Value',1);
    set(handles.startprofillastmat,'Value',0);


    %polynomials=equations('polynom',start,bvpfile,x1);
    %plotpoly(polynomials,x1,feval(bvpfile,'n'));
    valx1tau=equations('valx1tau',start,bvpfile,x1);
    x1tau=equations('x1tau',start,bvpfile,x1);
    n=feval(bvpfile,'n');

    ep=str2num(get(handles.edit24,'String'));


    figure
    
    if feval(bvpfile,'Infinite')&& feval(bvpfile,'EVP')==0


        if ep ~=0

            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end

            for i=1:n

                subplot(n,1,i);


                plot(tau_infinite,sol_infinite(i,:),'color','black');
                axis tight

                xlabel('t');
                ylabel(['z_',num2str(i)]);
            end

            title(subplot(n,1,1),'Initial values');
        else

            for i=1:n/2


                
                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),0);

            end

            for i=1:n/2

                subplot(n/2,1,i);


                plot(tau_infinite,sol_infinite(i,:),'color','black');
                axis tight
                xlabel('t');
                ylabel(['z_',num2str(i)]);

            end
            title(subplot(n/2,1,1),'Initial values');
        end

    elseif feval(bvpfile,'Infinite') && feval(bvpfile,'EVP')



        if ep == 0


            for i=1:n/2



                [sol_infinite(i,:),tau_infinite]= backtransf(valx1tau(i,2:end),valx1tau(i+n/2,2:end),x1tau(2:end),ep);



            end


            for i=1:n/2-2
%a=sol_infinite(n/2-1,1);
                subplot(n/2-2,1,i);
                plot(tau_infinite,sol_infinite(i,:),'color','black');
                axis tight
                xlabel('t');
                ylabel(['z_',num2str(i)]);



            end

title(subplot(n/2-2,1,1),'Initial values');

        else


            for i=1:n

                [sol_infinite(i,:),tau_infinite]= backtransf([],valx1tau(i,2:end),x1tau(2:end),ep);

            end


            for i=1:n-2
                subplot(n-2,1,i);



                plot(tau_infinite,sol_infinite(i,:),'color','black');
                
                a=sol_infinite(n-1,1);
axis tight
xlabel('t');
                ylabel(['z_',num2str(i)]);
            end

 title(subplot(n-2,1,1),'Initial values');
        end





    elseif feval(bvpfile,'Infinite') ==0 && feval(bvpfile,'EVP')



        for i=1:n-2
            subplot(n-2,1,i);
            plot(x1tau,valx1tau(i,:),'color','black');
            axis tight
xlabel('t');
                ylabel(['z_',num2str(i)]);
            a=valx1tau(n-1,1);
        end

title(subplot(n-2,1,1),'Initial values');
    elseif feval(bvpfile,'Infinite') ==0 && feval(bvpfile,'EVP')== 0
        
        for i=1:n
            subplot(n,1,i);
            plot(x1tau,valx1tau(i,:),'color','black');
            axis tight
xlabel('t');
                ylabel(['z_',num2str(i)]);
        end

title(subplot(n,1,1),'Initial values');
    end



    

  
 end

 

if get(handles.popupmenu1,'Value')==3 %3. use last solution

    datei=fopen('last.mat');
    if datei==-1
        err('bvps_errdlg24');
        return
    end
    fclose(datei);
    load last.mat
    parameter=str2num(get(handles.parameter,'String'));
    if parameter>0
        p=coeff(length(coeff)-parameter+1:length(coeff));
        p=p';
    end
    
    if  feval(bvpfile,'Infinite')
      help=strcat('[',num2str(tau_infinite,20),']');  
    else    
     help=strcat('[',num2str(x1tau,20),']');   
    end    
       
    set(handles.startprofilgitter,'String',help);
    
    
    
    if  feval(bvpfile,'Infinite')
        
       if feval(bvpfile,'EVP') 
           help='[';
           for i=1:size(sol_infinite,1)-2
               help=strcat(help,num2str(sol_infinite(i,:),20));
               help=strcat(help,';');
           end
           help=strcat(help,']');
       else
           
           help='[';
           for i=1:size(sol_infinite,1)
               help=strcat(help,num2str(sol_infinite(i,:),20));
               help=strcat(help,';');
           end
           help=strcat(help,']');
           
           
       end 
        
        

    else
        
        if feval(bvpfile,'EVP')

            help='[';
            for i=1:length(valx1tau(:,1))-2
                help=strcat(help,num2str(valx1tau(i,:),20));
                help=strcat(help,';');
            end
            help=strcat(help,']');

        else

            help='[';
            for i=1:length(valx1tau(:,1))
                help=strcat(help,num2str(valx1tau(i,:),20));
                help=strcat(help,';');
            end
            help=strcat(help,']');
        end
        
    end
    
    set(handles.startprofilwerte,'String',help);
    
    if  feval(bvpfile,'EVP')
      
 
      help=num2str(lambda);    
      set(handles.lambda,'String',help);
    
    end    
    
    
    
    if  feval(bvpfile,'Infinite')
      help=num2str(length(x1));    
    else    
       help=strcat('[',num2str(x1,10),']');    
    end    
   
    set(handles.neuesgitter,'String',help);
    if parameter>0
        help=strcat('[',num2str(p,20),']');
        set(handles.startprofilparameter,'String',help);
    else
        set(handles.startprofilparameter,'String','');
    end
    set(handles.checkbox2,'Value',1);
    set(handles.radiobutton4,'Value',0);
    set(handles.startprofillastmat,'Value',0);

    helpdlg('Modify "Alternative mesh" arbitrarily and push "SOLVE"!','Help');    

end


function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function popupmenu1_Callback(hObject, eventdata, handles)


function checkbox3_Callback(hObject, eventdata, handles)


function checkbox4_Callback(hObject, eventdata, handles)
if get(handles.checkbox4,'Value')==1
    helpdlg('Fill in the fields "Initial mesh" and "Initial values / parameters", they will be saved in the bvpfile. The field "Alternative mesh" is optional, because the values from the field "mesh" will be used when left out. If no inital mesh is inserted the inital profile is set to the constant 1!','Help');
    set(handles.edit2,'String','1');
end


function graphischdarstellen_Callback(hObject, eventdata, handles)
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


if get(mainGUIdata.radiobutton12,'Value')==1
    
    %feval(bvpfile,'Infinite')
    
    cd(aktuellesverzeichnis);
    plotrange
    cd(arbeitsverzeichnis)  


 

elseif feval(bvpfile,'EVP') && get(mainGUIdata.radiobutton12,'Value')==0
   
    cd(aktuellesverzeichnis);
    plot_results(x1tau,valx1tau,n-2,1,0);



else
   
    
    cd(aktuellesverzeichnis);
    plot_results(x1tau,valx1tau,n,0,0);

     
end




cd(aktuellesverzeichnis);



% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


function pushbutton29_Callback(hObject, eventdata, handles)
load last.mat
bvpfile=get(handles.edit1,'String');
pfad=get(handles.pfad,'String');
bvpfilem=bvpfile;
bvpfile=strrep(bvpfile,'.m','');
datei=fopen(strcat(pfad,bvpfilem),'r');
if datei==-1
    err('bvps_errdlg03');
    return;
end
inhaltdateibak=fread(datei);
fclose(datei);
%zufallsname=strcat('temp_bvpsuite_file_',strcat(strcat(strrep(strrep(num2str(clock),' ',''),'.',''),strrep(num2str(randsrc(1,5,[0 1 2 3 4 5 6 7 8 9])),' ','')),'.m'));
%zufallsname='temp_bvpsuite_file.m';
%delete('temp_bvpsuite_file_*.m');
%datei=fopen(zufallsname,'w');
%fprintf(datei,'%s',char(inhaltdateibak));
%fclose(datei);
%zufallsname2=strrep(zufallsname,'.m','');
%Dummyauswertung
feval(bvpfile,'x1');
%here you get the vector of points from the GUI
t=str2num(get(handles.edit21,'String'));
if length(t)==0
    errordlg('Fill in the field!','Error');
    return
end



if feval(bvpfile,'Infinite') && feval(bvpfile,'EVP')

    for k=1:length(t)

        if (t(k)< x1(1))
            help1=num2str(x1(1));
            intervall=strcat(' [',help1,',','infinity',')');
            text=strcat('The value is not in the calculated range!',intervall,'!');
            errordlg(text,'Error');
            return
        end

        if feval(bvpfile,'Endpoint')==0

            if t(k)>=1
                value(:,k)=equations('wert',coeff,bvpfile,x1,1/t(k));
                wert(:,k)=value(feval(bvpfile,'n')/2+1:end-2,k);

            else
                value(:,k)=equations('wert',coeff,bvpfile,x1,t(k));
                wert(:,k)=value(1:feval(bvpfile,'n')/2-2,k);

            end
        else

            value(:,k)=equations('wert',coeff,bvpfile,x1,1/t(k));
            wert(:,k)=value(1:end-2,k);



        end

    end





    assignin('base','t',t);
    assignin('base','wert',wert);




elseif feval(bvpfile,'Infinite') && feval(bvpfile,'EVP')==0

    for k=1:length(t)

        if (t(k)< x1(1))
            help1=num2str(x1(1));
            intervall=strcat(' [',help1,',','infinity',')');
            text=strcat('The value is not in the calculated range!',intervall,'!');
            errordlg(text,'Error');
            return
        end

        
        if feval(bvpfile,'Endpoint')==0
        
         if t(k)>=1
               value(:,k)=equations('wert',coeff,bvpfile,x1,1/t(k));
              wert(:,k)=value(feval(bvpfile,'n')/2+1:end,k);

         else
              value(:,k)=equations('wert',coeff,bvpfile,x1,t(k));
             wert(:,k)=value(1:feval(bvpfile,'n')/2,k);

         end
         
        else
            
             value(:,k)=equations('wert',coeff,bvpfile,x1,t(k));
             wert(:,k)=value(:,k); 
            
            
        end  
         
    end




    assignin('base','t',t);
    assignin('base','wert',wert);




elseif feval(bvpfile,'EVP') && feval(bvpfile,'Infinite')==0




    for k=1:length(t)

        if (t(k)< x1(1) || t(k)> x1(end))
            help1=num2str(x1(1));
            help2=num2str(x1(end));
            intervall=strcat(' [',help1,',',help2,']');
            text=strcat('The value is not in the calculated range!',intervall,'!');
            errordlg(text,'Error');
            return
        end

        value(:,k)=equations('wert',coeff,bvpfile,x1,t(k));
        wert(:,k)=value(1:(feval(bvpfile,'n')-2),k);

    end


    assignin('base','t',t);
    assignin('base','wert',wert);

else

    for k=1:length(t)

        if (t(k)< x1(1) || t(k)> x1(end))
            help1=num2str(x1(1));
            help2=num2str(x1(end));
            intervall=strcat(' [',help1,',',help2,']');
            text=strcat('The value is not in the calculated range!',intervall,'!');
            errordlg(text,'Error');
            return
        end
    end

    wert =equations('wert',coeff,bvpfile,x1,t);
    assignin('base','t',t);
    assignin('base','wert',wert);

end

open t;
open wert;
%help='[';
%for i=1:length(wert(:,1))
%    help=strcat(help,num2str(wert(i,:)));
%    help=strcat(help,';');
%end
%help=strcat(help,']');
%set(handles.edit22,'String',help);



function edit22_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit22_Callback(hObject, eventdata, handles)

function pushbutton31_Callback(hObject, eventdata, handles)
close



% --- Executes on button press in einstellungen.
function einstellungen_Callback(hObject, eventdata, handles)
settings









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



function pfad_Callback(hObject, eventdata, handles)
% hObject    handle to pfad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pfad as text
%        str2double(get(hObject,'String')) returns contents of pfad as a double


% --- Executes during object creation, after setting all properties.
function parameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function parameter_Callback(hObject, eventdata, handles)
% hObject    handle to parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parameter as text
%        str2double(get(hObject,'String')) returns contents of parameter as a double


% --- Executes during object creation, after setting all properties.
function startprofilparameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startprofilparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function startprofilparameter_Callback(hObject, eventdata, handles)
% hObject    handle to startprofilparameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startprofilparameter as text
%        str2double(get(hObject,'String')) returns contents of startprofilparameter as a double


% --- Executes on button press in konvergenzord.





% --- Executes on button press in delaltgitter.
function delaltgitter_Callback(hObject, eventdata, handles)
% hObject    handle to delaltgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.neuesgitter,'String','');


% --- Executes on button press in delstartprofilwerte.
function delstartprofilwerte_Callback(hObject, eventdata, handles)
% hObject    handle to delstartprofilwerte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.startprofilwerte,'String','');
set(handles.startprofilparameter,'String','');
set(handles.lambda,'String','');


% --- Executes on button press in delstartprofilgitter.
function delstartprofilgitter_Callback(hObject, eventdata, handles)
% hObject    handle to delstartprofilgitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.startprofilgitter,'String','');




function variablen_Callback(hObject, eventdata, handles)
% hObject    handle to variablen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of variablen as text
%        str2double(get(hObject,'String')) returns contents of variablen as a double


% --- Executes during object creation, after setting all properties.
function variablen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to variablen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function c_Callback(hObject, eventdata, handles)
% hObject    handle to c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c as text
%        str2double(get(hObject,'String')) returns contents of c as a double


% --- Executes during object creation, after setting all properties.
function c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in fehler.
function fehler_Callback(hObject, eventdata, handles)
bvpfile=get(handles.edit1,'String');
bvpfilem=bvpfile;
bvpfile=strrep(bvpfile,'.m','');
pfad=get(handles.pfad,'String');
%zufallsname=strcat('temp_bvpsuite_file_',strcat(strcat(strrep(strrep(num2str(clock),' ',''),'.',''),strrep(num2str(randsrc(1,5,[0 1 2 3 4 5 6 7 8 9])),' ','')),'.m'));
%zufallsname='temp_bvpsuite_file.m';
datei=fopen(strcat(pfad,bvpfilem),'r');
if datei==-1
    errordlg('The file does not exist. Check its name.','Error');
    return;
end
inhaltdatei=fread(datei);
fclose(datei);
inhaltdateibak=inhaltdatei;
inhaltdatei=char(inhaltdatei');
if (length(strfind(inhaltdatei,'Kontrollnummer43753976430976655143'))>0)
    err('bvps_errdlg28');
    return;
end
gueltig=strfind(inhaltdatei,'Kontrollnummer43753976430976655144');
if (length(gueltig)>0) 
      err('bvps_errdlg29');
      return;
end

if (length(gueltig)==0 & length(strfind(inhaltdatei,'Kontrollnummer43753976430976655145'))==0  && length(strfind(inhaltdatei,'Kontrollnummer43753976430976655143'))==0)
      errordlg('The chosen file is not compatible with any version of bvpsuite!','Error');
      return;
end
checkm3=strfind(bvpfilem,'bvpsuite.m');
if(length(checkm3)~=0)
    errordlg('bvpsuite is not a bvpfile!','Error');
    return;
end
delete('temp_bvpsuite_file_*.m');
%datei=fopen(zufallsname,'w');
%fprintf(datei,'%s',char(inhaltdateibak));
%fclose(datei);
%zufallsname2=strrep(zufallsname,'.m','');


load last.mat
load options.mat
[x1tau,valerror,maxfehlerwerte,koeff_2,x1tau_2,wert_2x1undtau_2,x1_2,polynome_2,wert2_x1_2,valerror2,tau_infinite_2,error_infinite,error_infinite_2,sol_infinite_2,lambda_2,eigenfunction_2]=errorestimate(coeff,bvpfile,1,x1,bvpopt,ausgabe);

% open x1tau;
% open valerror;


n=feval(bvpfile,'n');

if feval(bvpfile,'Infinite') 
    
    
    if feval(bvpfile,'EVP')
           
        
        assignin('base','tau_infinite',tau_infinite);
        assignin('base','error_infinite',error_infinite);
        assignin('base','tau_infinite_2',tau_infinite_2);
        assignin('base','error_infinite_2',error_infinite_2);
        coeff=koeff_2;
        x1=x1_2;
        valx1=wert2_x1_2;
        x1tau=x1tau_2;
        valx1tau=wert_2x1undtau_2;
        polynomials=polynome_2;
        tau_infinite=tau_infinite_2;
        sol_infinite=sol_infinite_2;
        lambda=lambda_2;
        eigenfunction=eigenfunction_2;
        save lasterr coeff x1 valx1 x1tau valx1tau polynomials tau_infinite sol_infinite lambda eigenfunction;
    else
        
        assignin('base','tau_infinite',tau_infinite);
        assignin('base','error_infinite',error_infinite);
        assignin('base','tau_infinite_2',tau_infinite_2);
        assignin('base','error_infinite_2',error_infinite_2);
        coeff=koeff_2;
        x1=x1_2;
        valx1=wert2_x1_2;
        x1tau=x1tau_2;
        valx1tau=wert_2x1undtau_2;
        polynomials=polynome_2;
        tau_infinite=tau_infinite_2;
        sol_infinite=sol_infinite_2;
        lambda=lambda_2;
        eigenfunction=eigenfunction_2;
        save lasterr coeff x1 valx1 x1tau valx1tau polynomials tau_infinite sol_infinite;
    end 
    
else
    
    if feval(bvpfile,'EVP')
        assignin('base','x1tau',x1tau);
        assignin('base','valerror',valerror(1:n-2,:));
        assignin('base','x1tau_2',x1tau_2);
        assignin('base','valerror2',valerror2(1:n-2,:));
        coeff=koeff_2;
        x1=x1_2;
        valx1=wert2_x1_2;
        x1tau=x1tau_2;
        valx1tau=wert_2x1undtau_2;
        polynomials=polynome_2;
        tau_infinite=tau_infinite_2;
        sol_infinite=sol_infinite_2;
        lambda=lambda_2;
        eigenfunction=eigenfunction_2;
        
        
        save lasterr coeff x1 valx1 x1tau valx1tau polynomials lambda sol_infinite eigenfunction;

    else

        assignin('base','x1tau',x1tau);
        assignin('base','valerror',valerror);
        assignin('base','x1tau_2',x1tau_2);
        assignin('base','valerror2',valerror2);
        coeff=koeff_2;
        x1=x1_2;
        valx1=wert2_x1_2;
        x1tau=x1tau_2;
        valx1tau=wert_2x1undtau_2;
        polynomials=polynome_2;
        tau_infinite=tau_infinite_2;
        sol_infinite=sol_infinite_2;
        lambda=lambda_2;
        eigenfunction=eigenfunction_2;
        
        
        save lasterr coeff x1 valx1 x1tau valx1tau polynomials;
        
    end 
    
    
end 






function meshadaptation_Callback(hObject, eventdata, handles)

%  if get(handles.meshadaptation,'Value')==1
%       set(handles.meshadaptation2,'Value',0);
%  end
 






function startprofildatei_Callback(hObject, eventdata, handles)
% hObject    handle to startprofildatei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startprofildatei as text
%        str2double(get(hObject,'String')) returns contents of startprofildatei as a double


% --- Executes during object creation, after setting all properties.
function startprofildatei_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startprofildatei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in pathfollowing.
function pathfollowing_Callback(hObject, eventdata, handles)


%if get(handles.radiobutton12,'Value')

%helpdlg('Your problem is posed on a semi-infinite interval. See the bvpsuite manual for instructions how to use the pathfollowing option for this problem class.','help');
    
%else
berechnen(4,hObject,eventdata,handles);
%end 

%bvpsuite question mark help-texts

function help1_Callback(hObject, eventdata, handles)
helpdlg('Choose the path and input the name of the m-file (input: *.m), which you would like to edit or last.log if you would like to see your last input. Use "File --> Open example" to find your file.','Help');

function help2_Callback(hObject, eventdata, handles)
helpdlg('The program solves the collocation-equation-system by Newton''s iteration. It is therefore necessary, to chose an initial value. All coefficients will be set to this decimal value, which can only be proposed for examples, that are not very sensitive to these values. Fill in the fields "Initial Grid", "Initial Values" (respectivly "Parameters") and chose "Save initial profile", to achieve better convergence. The value in the field ''Initial Value'' will be ignored in your next calculation after pushing the button ''Save''.','Help');

function help3_Callback(hObject, eventdata, handles)
helpdlg('Choose the orders of the solution components in the form of a Matlab row vector. (Be careful to choose the correct order!) E.g.: [3 4 2], if you have 3 equations and therefore 3 solution components with the highest derivative(z1)=3, derivative(z2)=4 and derivative(z3)=2, respectively. Specify the number of unknown parameters, denote them by p1,p2,..., in the equations. When solving a BVP posed on the interval [a, b] with two-point boundary conditions specified at a and b, leave the field "c" empty. Note that bvpsuite is also capable of solving systems where the boundary conditions are posed in the interior of the interval [a, b]. If you solve such a problem, specify a Matlab row vector with the positions of the interior points. In the field "Boundary / Additional conditions" specify z2(c1)=...,z1(c3)=..., instead of z2(a)=.... The variables a and b for the endpoints of the interval may no longer be used!','Help')

function help4_Callback(hObject, eventdata, handles)
helpdlg('The number of predefined collocation points (Gaussian, Lobatto, or Uniform) has to be specified in "Number". When the button "User" is ticked, the user can specify the number and his own distribution of collocation points. When the number of collocation points is m, then one has to describe their distribution in the field "Partition rho_i". There, m values 0 <= rhi_i <= 1, i = 1, . . . ,m in form of a Matlab row vector have to be specified. They describe the positions of the collocation points in the interval [0, 1].','Help');

function help5_Callback(hObject, eventdata, handles)
helpdlg('Specify the number of predifened collocation-points (Gauss, Lobatto, equidistant) or chose your own Matlab-row-vector of collocation-points when radio button ''User'' is checked.','Help');

function help6_Callback(hObject, eventdata, handles)
helpdlg('Specify the position of the mesh points, e.g. [0 1.7 2.5 3]. Note that a and b, the left and right endpoints of the interval have to be in the list. In this case a=0 and b=3. You can also specify the mesh in the form a:h:b, where h is the stepsize. In case that your problem is posed on a semi-infinite interval only insert the number of subintervals, e.g. 50. This defines an equidistant mesh on the finite domain. Later, a graphical back transformation of the solution is carried out by the code.','Help')

function help7_Callback(hObject, eventdata, handles)
helpdlg('You can define string replacements for parameters frequently occurring in the equations and boundary conditions. Syntax: string_to_be_replaced_1 = replacing_string_1; string_to_be_replaced_2 = replacing_string_2; ... E.g.: m = 5*z1''''-8; j=3; Be careful: The strings t,z,z1,z2,...,p1,p2,..., and a, b must not be replaced. Be careful with string replacements var=4; var1=7; because the second string will be replaced by 41=7.','Help');

function help8_Callback(hObject, eventdata, handles)
helpdlg('Specify the differential equations using variables z_1(t), ..., z_n(t), z''1(t), ..., z''_n(t), z''''_1 (t), ..., z''''_n(t), written as z1,...,zn,z1'',...,zn'',z1'''',...,zn''''. Higher derivatives (>=3), e.g. z(3), have to be written as z1d3. The equations can be specified in the form left side 1=right side 1; left side 2=right side 2. E.g.: z1''''=z2;z2''''=-z1. If you want to solve an eigenvalue problem always denote the eigenvalue by lambda.','Help');

function help9_Callback(hObject, eventdata, handles)
helpdlg('Specify the boundary / additional conditions separated by semicolons. The number of necessary conditions is equal to the sum of orders and parameters. Write z_i(a)=: zi(a), z''_i(a) =: zi''(a), z_i(b) =: zi(b), z''_i(b) =: zi''(b). For derivatives k >= 3 type zidk(a). E.g.: z1(a) = 3; z2(a) = 0; z1(b) = 0; z3d4(b)=7. Specify the boundary conditions at the interior points in the following way: z2(c1)=...;z1(c3)=... Attention: Here, a and b must not be used!','Help');

function help10_Callback(hObject, eventdata, handles)
helpdlg('Specify the mesh points for the initial solution values, in form of a Matlab row vector, e.g. [3 4.5 7].','Help');

function help11_Callback(hObject, eventdata, handles)
helpdlg('Choose the values of the initial profile in the form of a Matlab matrix. The rows in this matrix contain the guess for the solution values at the mesh points, e.g. for 3 mesh points and two solution components, [1 4 2;3 6 8]. If you leave the field "Initial values" empty, the initial profile will be set to the constant function equal to 1. In this case "Initial mesh" also has to be left empty. Additionally, for unknown parameters you can also specify their initial values in the corresponding field. For the solution of EVPs you can provide a guess for the eigenvalue ("lambda") and for the eigenfunction ("Initial values"). Alternatively, you can insert a value for "lambda" only. Then the initial profile for the eigenfunction is set to the constant function equal to 1.','Help');

function help12_Callback(hObject, eventdata, handles)
helpdlg('This mesh replaces the "Mesh" on the left hand side of the GUI window if the option "Profile in the fields above" or the option "As profile" is ticked.','Help');

function help13_Callback(hObject, eventdata, handles)
helpdlg('Executes the selected action in the drop down menu. "Plot initial profile above" shows a graphical output of the profile that was defined in the fields above. "Plot saved initial profile" shows the initial profile which was saved together with the file. "Use last solution as initial profile" replaces the values above by the solution computed in the last run.','Help');

function help14_Callback(hObject, eventdata, handles)
helpdlg('"Saved Profile": Solve the problem using the saved initial profile. "Profile in the fields above": Use the initial profile displayed in the fields above for the next run. "As profile": Use solution values stored in "last" (or in another *.mat file from the current directory) on the alternative mesh as initial guess for the next run. Attention: Do not change the number of collocation points when using the third option!','Help')

function help15_Callback(hObject, eventdata, handles)
helpdlg('Starts the computations. Mesh adaptation is the default computational mode and therefore the check box "Meshadaptation" is checked when the GUI window opens. The mesh adaptation algorithm is based on the residual and global error control.','Help');

function help16_Callback(hObject, eventdata, handles)
helpdlg('The results of the last calculation (stored in "last.mat") are exported to the workspace (overwriting existing variables with the same names). "coeff": The coefficients with respect to the Runge-Kutta basis for the collocation polynomials and numerical values of the unknowns parameters. "x1": The mesh points. "x1tau": The mesh points and the collocation points. "valx1tau": The values of the numerical solution in the mesh and collocation points. "parameters": Numerical approximation for the parameters. In the context of problems posed on semi-infinite interval, "tau_infite" is "x1tau" transformed to a truncated interval. For an EVP posed on a semi-infinite domain the "eigenfunction" is an approximation of the eigenfunction on "tau_infite", or otherwise, on "x1tau". For a BVP posed on a semi-infinite interval "sol_infinite" is an approximation of its solution on "tau_infinite".','Help');

function help17_Callback(hObject, eventdata, handles)
helpdlg('Plots the results stored in "last.mat" on "x1tau" or "tau_infinite".','Help');

function help18_Callback(hObject, eventdata, handles)
helpdlg('Input for a Matlab row vector containing additional, user selected points, for the values of the numerical approximation at those points.','Help');





% --- Executes on button press in help20.
function help20_Callback(hObject, eventdata, handles)
% hObject    handle to help20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Tick this check box to declare that your example constitutes an eigenvalue problem.','Help');

% --- Executes on button press in help21.
function help21_Callback(hObject, eventdata, handles)
% hObject    handle to help21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Tick this check box in case that your problem is posed on a semi-infinite interval [a,inifinity), a >= 0 and type the left endpoint of the interval in the edit field to the right. A transformation of the problem is carried out in order to transform it to a finite interval. As a result, the code computes a solution on this finite domain. A graphical back transformation allocates the solution on a truncated interval [a,L]. The right point L has to be specified by the user, cf. "Plot".','Help');

% Hint: get(hObject,'Value') returns toggle state of radiobutton15
% --- Executes on button press in help22.
function help22_Callback(hObject, eventdata, handles)
% hObject    handle to help22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


helpdlg(strcat('Check this button if your problem is of the form:',strcat('-(p*z''',')''','+q*y= \lambda*g*y'),strcat('. Insert p, q, and g in the following fields.')),'Help');

% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Enter the left endpoint of the semi-infinite interval. E.g. 0.','Help');




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


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --- Executes on button press in help20.
%function help20_Callback(hObject, eventdata, handles)
% hObject    handle to help20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press over help4 with no controls selected.
function help4_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to help4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over help20.
function help20_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to help20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over help3.
function help3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to help3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton12














function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Check this button in case that you want to calculate several eigenvalues (EV) and eigenfunctions (EF) of your problem.','Help');










% --- Executes on button press in meshadaptation2.
function meshadaptation2_Callback(hObject, eventdata, handles)
% hObject    handle to meshadaptation2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of meshadaptation2




% --- Executes on button press in automatic.
function automatic_Callback(hObject, eventdata, handles)
% hObject    handle to automatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of automatic

%Info: automatic order is not changed if no filename is given!!!!

if get(handles.automatic,'Value')==0
   set(handles.edit13,'BackgroundColor','white');
   set(handles.edit13,'Enable','on');
   set(handles.edit13,'String','');
   
    helpdlg('Press "Save" to save changes!','help');

elseif get(handles.automatic,'Value')==1
  
    set(handles.edit13,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   set(handles.edit13,'Enable','off');    

  load options.mat 
  
  %open abstol
  
   
  if  min(abstolgitter) <= 10^-6                             

    set(handles.edit13,'String',8);
  
  elseif  min(abstolgitter) > 10^-6 &&  min(abstolgitter) <= 10^-4
      
    set(handles.edit13,'String',6);
  
  elseif  min(abstolgitter) > 10^-4 &&  min(abstolgitter) <= 10^-2
      
    set(handles.edit13,'String',4);  
 
  elseif  min(abstolgitter) > 10^-2 
  
    set(handles.edit13,'String',2);     
  end     
    
    
  helpdlg('Press "Save" to save changes!','help');

  
end



% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


helpdlg('When the GUI window opens this check box is ticked. In such a case, the code automatically determines the number of collocation points to be used for the computations depending on the prescribed "Absolute tolerance" in "Settings". Note that you have to press "Save" again after you have ticked "Automatic" in order to save the changes for the computation. For help on "Number / Partition rho_i" refer to the previous "?".','Help');



function lambda_Callback(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda as text
%        str2double(get(hObject,'String')) returns contents of lambda as a double


% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aktuellesverzeichnis=cd;
arbeitsverzeichnis=get(handles.pfad,'String');
cd(arbeitsverzeichnis);

[FileName,PathName,index] = uigetfile('*.m;*.log','Choose the m-file');

if index>0
    set(handles.edit1,'String',FileName);
    set(handles.pfad,'String',PathName);
end
cd(aktuellesverzeichnis);
if index>0
    oeffnen_Callback(hObject, eventdata, handles);
end




function radiobutton1_Callback(hObject, eventdata, handles)
if get(handles.radiobutton1,'Value')==0
    set(handles.radiobutton1,'Value',1);
else
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);
    set(handles.lobatto,'Value',0);
end




% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Matlab Code BVPSUITE1.1 (Release July 2010). Authors: Georg Kitzhofer, Othmar Koch, Gernot Pulverer, Christa Simon and Ewa B. Weinmueller. Institute for Analysis and Scientific Computing, Vienna University of Technology, Vienna, Austria. BVPSUITE aims for the efficient numerical solution of boundary value problems in ordinary differential equations. In its scope are fully implicit problems of mixed orders, parameter dependent problems, problems with unknown parameters, problems posed on semi-infinite intervals, eigenvalue problems and differential algebraic equations of index-1. Singularities of the first and second kind as well as the space singularities arising in the differential operator are admissible. The singular points may occur at either one endpoint or at both endpoints of the interval of integration. The program can be applied directly to such singular boundary value problems and no pre-handling is necessary. Boundary conditions can be also posed in the interior points of the interval. For higher efficiency an estimate of the global error and an adaptive mesh selection is provided. Automatic transformation of problems posed on semi-infinite intervals to a finite domain makes the solution of such problems also accessible to the code. For more information and for the manual contact http://www.math.tuwien.ac.at/~ewa or e.weinmueller@tuwien.ac.at','Info');

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on key press over oeffnen with no controls selected.
function oeffnen_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to oeffnen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('Plots a graphical output of the global discretization error and copies the values to the Matlab workspace. The out- put variables are "x1tau", "x1tau2", and "valerror", "valerror2". There are two error plots. In the Figure 1, the values in "x1tau" and "valerror" are associated with the coarse mesh and show the error estimate of the absolute error of the solution computed on the coarse mesh. In the Figure 2, the values "x1tau2" and "valerror2" are associated with the fine mesh (with twice as many points) and show the estimate of the absolute error of the solution computed on the fine mesh. These two meshes are used during the error estimation proce- dure. In case of a problem posed on the semi-infinite interval the output variables are "tau_infinite", "error_infinite", and "tau_infinite_2", "error_infinite_2".','help');




function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
