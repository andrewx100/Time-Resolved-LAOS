function varargout = postAnal(varargin)
% POSTANAL MATLAB code for postAnal.fig
%      POSTANAL, by itself, creates a new POSTANAL or raises the existing
%      singleton*.
%
%      H = POSTANAL returns the handle to a new POSTANAL or the handle to
%      the existing singleton*.
%
%      POSTANAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTANAL.M with the given input arguments.
%
%      POSTANAL('Property','Value',...) creates a new POSTANAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before postAnal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to postAnal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help postAnal

% Last Modified by GUIDE v2.5 03-Sep-2011 23:12:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @postAnal_OpeningFcn, ...
    'gui_OutputFcn',  @postAnal_OutputFcn, ...
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


% --- Executes just before postAnal is made visible.
function postAnal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to postAnal (see VARARGIN)

% Choose default command line output for postAnal
handles.output = hObject;
handles.blnFileSelected = false;
handles.blnValidFreq = false;
handles.strLastPath='\';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes postAnal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = postAnal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txtFileName_Callback(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFileName as text
%        str2double(get(hObject,'String')) returns contents of txtFileName as a double


% --- Executes during object creation, after setting all properties.
function txtFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnBrowse.
function btnBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName] = uigetfile('*.txt','Select a data file to process',handles.strLastPath);
if ischar(FileName)
    handles.strLastPath = PathName;
    handles.strFileName = [PathName FileName];
    set(handles.txtFileName, 'String',[PathName FileName]);
    % handles.x=uiimport(handles.strFileName);
    fid=fopen(handles.strFileName);
    dataHead = fscanf(fid, '%d %d %d',3);
    handles.pts_per_channel=dataHead(1);
    handles.nr_oversampl=dataHead(2);
    handles.freq_max=dataHead(3);
    handles.blnFileSelected=true;
    guidata(hObject,handles);
else
    msgbox ('No file selected.','Notice','warn','modal')
    set(handles.txtFileName,'String','Enter a file name.');
    handles.blnFileSelected=false;
    set(handles.btnProcess,'Enable','off');
    guidata(hObject, handles);
end
if handles.blnFileSelected && handles.blnValidFreq
    set(handles.btnProcess,'Enable','on');
    guidata(hObject, handles);
else
    set(handles.btnProcess,'Enable','off');
    guidata(hObject, handles);
end
guidata(hObject, handles);



% --- Executes on button press in btnProcess.
function btnProcess_Callback(hObject, eventdata, handles)
% hObject    handle to btnProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.strFileName
A=importdata(handles.strFileName);
x=A.data;
[row,col]=size(x);

% Signal and sampling parameters
fi = handles.dblFI;
fs_prime=handles.pts_per_channel/handles.nr_oversampl;


% Add a time column
t=0:1/fs_prime:(row-1)/fs_prime;
t=t';
x=[t x];

% Total number of cycles
Nr_of_Cycle=int32(t(row)/fi);

% Rows per cycle
row_per_cycle = round(fs_prime/fi);
data=[];

% ti=0;
% data_1=x(1:3*row_per_cycle,:);
% [row_1,col_1]=size(data_1);
%
% freq_bin=0:row_1-1;
% freq_res=fs_prime/row_1;
% freq_1=freq_bin * freq_res;
% Xraw=fft(data_1(:,3))/row_1;
% cutoff=ceil(row_1/2);
% freq_1=freq_1(1:cutoff);
% Xraw=Xraw(1:cutoff);
% % plot(freq_1,abs(Xraw))
%     f1 = round(fi / freq_res);
%     f3 = round(3 * fi / freq_res);
%     f5 = round(5 * fi / freq_res);
%     f7 = round(7 * fi / freq_res);
%     f9 = round(9 * fi / freq_res);
%         I31 = abs(Xraw(f3+1))/abs(Xraw(f1+1));
%     I51 = abs(Xraw(f5+1))/abs(Xraw(f1+1));
%     I71 = abs(Xraw(f7+1))/abs(Xraw(f1+1));
%     I91 = abs(Xraw(f9+1))/abs(Xraw(f1+1));
%
% data = [ti I31 I51 I71 I91]


for i = 0:1:Nr_of_Cycle-3;
    begin_row = round(i * fs_prime / fi)+1;
    ti=x(begin_row,1);
    data_i = x(begin_row:begin_row + 2 * row_per_cycle,:); % Take exactly 2 cycles
    %     plot(data_i(:,1),data_i(:,2))
    [row_i, col_i]=size(data_i);
    
    freq_bin = 0:row_i-1; % vector of frequency bins
    freq_res = fs_prime/row_i; % frequency resolution
    freq_i = freq_bin * freq_res; % frequency axis
    Xraw = fft(data_i(:,3))/row_i; % normalized fft
    cutoff = ceil(row_i/2); % only use the first half of the FFT
    freq_i=freq_i(1:cutoff);
    Xraw=Xraw(1:cutoff);
    %     plot(freq_i,abs(Xraw))
    % frequency bins of the harmonics
    f1 = round(fi / freq_res);
    f3 = round(3 * fi / freq_res);
    f5 = round(5 * fi / freq_res);
    f7 = round(7 * fi / freq_res);
    f9 = round(9 * fi / freq_res);
    
    % Relative intensity of harmonics
    I31 = abs(Xraw(f3+1))/abs(Xraw(f1+1));
    I51 = abs(Xraw(f5+1))/abs(Xraw(f1+1));
    I71 = abs(Xraw(f7+1))/abs(Xraw(f1+1));
    I91 = abs(Xraw(f9+1))/abs(Xraw(f1+1));
    
    % MITlaos
    [G_M,G_L,Eta_M,Eta_L,NE,NV]=MITLaos(fi,fs_prime,300,data_i(:,2),data_i(:,3));
    
    if i==0
        data=[ti I31 I51 I71 I91 G_M G_L Eta_M Eta_L NE NV];
    else
        data=[data;ti I31 I51 I71 I91 G_M G_L Eta_M Eta_L NE NV];
    end
  
    
    
end
newFileName=get(handles.txtFileName, 'String');
newFileName=[newFileName(1:length(newFileName)-4) 'ft.txt'];
size(data);
save(newFileName,'data','-ascii');
set(handles.txtFileName,'String','Enter a file name.');
set(handles.txtFI,'String','0');
handles.blnValidFreq=false;
handles.blnFileSelected=false;
set(handles.btnProcess,'Enable','off');
guidata(hObject, handles);



function txtFI_Callback(hObject, eventdata, handles)
% hObject    handle to txtFI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFI as text
%        str2double(get(hObject,'String')) returns contents of txtFI as a double
handles.dblFI = str2double(get(hObject,'string'));

% Check the validity of the entered content
if isnan(handles.dblFI) || (handles.dblFI <= 0)
    errordlg('You must enter a positive, numeric value','Bad Input','modal')
    set(hObject, 'String', '1');
    set(handles.btnProcess,'Enable','off');
    handles.blnValidFreq=false;
    guidata(hObject, handles);
    uicontrol(hObject)
    return
end
handles.blnValidFreq=true;
handles.blnFileSelected
handles.blnFileSelected & handles.blnValidFreq
if handles.blnFileSelected && handles.blnValidFreq
    set(handles.btnProcess,'Enable','on');
    guidata(hObject, handles);
else
    set(handles.btnProcess,'Enable','off');
    guidata(hObject, handles);
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function txtFI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
