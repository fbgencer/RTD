function varargout = gui101ders(varargin)
% GUI101DERS MATLAB code for gui101ders.fig
%      GUI101DERS, by itself, creates a new GUI101DERS or raises the existing
%      singleton*.
%
%      H = GUI101DERS returns the handle to a new GUI101DERS or the handle to
%      the existing singleton*.
%
%      GUI101DERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI101DERS.M with the given input arguments.
%
%      GUI101DERS('Property','Value',...) creates a new GUI101DERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui101ders_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui101ders_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui101ders

% Last Modified by GUIDE v2.5 13-May-2018 16:16:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui101ders_OpeningFcn, ...
                   'gui_OutputFcn',  @gui101ders_OutputFcn, ...
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


% --- Executes just before gui101ders is made visible.
function gui101ders_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui101ders (see VARARGIN)
handles.peaks=peaks(35);
handles.membrane=membrane;

[x,y]=meshgrid(-8:0.5:8);
r=sqrt(x.^2+y.^2)+eps;
sinc=sin(r)./r;
handles.sinc=sinc;
handles.current_data=handles.sinc;
%surf(handles.current_data);


% Choose default command line output for gui101ders
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui101ders wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui101ders_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in peakbutton.
function peakbutton_Callback(hObject, eventdata, handles)
% hObject    handle to peakbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns peakbutton contents as cell array
%        contents{get(hObject,'Value')} returns selected item from peakbutton
deger=get(hObject,'Value');
satir=get(hObject,'String');
switch satir{deger}
    case 'Data1'
        handles.current_data=handles.peaks;
    case 'Data2'
        handles.current_data=handles.membrane;
    case 'Data3'
        handles.current_data=handles.sinc;  
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function peakbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','red');
end


% --- Executes on button press in surfbutton.
function surfbutton_Callback(hObject, eventdata, handles)
% hObject    handle to surfbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
surf(handles.current_data);

% --- Executes on button press in meshbutton.
function meshbutton_Callback(hObject, eventdata, handles)
% hObject    handle to meshbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mesh(handles.current_data);

% --- Executes on button press in counterbutton.
function counterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to counterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contour(handles.current_data);


% --- Executes on button press in klik.
function klik_Callback(hObject, eventdata, handles)
% hObject    handle to klik (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of klik
tik=get(hObject,'Value');
switch tik
    case 1
        view(2);
    case 0
        view(3);
end;


% --- Executes during object creation, after setting all properties.
function Resim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Resim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Resim
