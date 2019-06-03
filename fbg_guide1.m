function varargout = fbg_guide1(varargin)
% FBG_GUIDE1 MATLAB code for fbg_guide1.fig
%      FBG_GUIDE1, by itself, creates a new FBG_GUIDE1 or raises the existing
%      singleton*.
%
%      H = FBG_GUIDE1 returns the handle to a new FBG_GUIDE1 or the handle to
%      the existing singleton*.
%
%      FBG_GUIDE1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FBG_GUIDE1.M with the given input arguments.
%
%      FBG_GUIDE1('Property','Value',...) creates a new FBG_GUIDE1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fbg_guide1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fbg_guide1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fbg_guide1

% Last Modified by GUIDE v2.5 13-May-2018 16:32:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fbg_guide1_OpeningFcn, ...
                   'gui_OutputFcn',  @fbg_guide1_OutputFcn, ...
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


% --- Executes just before fbg_guide1 is made visible.
function fbg_guide1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fbg_guide1 (see VARARGIN)
handles.peaks=peaks(35);
handles.membrane=membrane;
[x,y] = meshgrid(-8:.5:8);
r = sqrt(x.^2+y.^2) + eps;
sinc = sin(r)./r;
handles.sinc = sinc;
% Set the current data value.
handles.current_data = handles.peaks;
surf(handles.current_data);
% Choose default command line output for fbg_guide1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fbg_guide1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Create the data to plot.


% --- Outputs from this function are returned to the command line.
function varargout = fbg_guide1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in surf_button.
function surf_button_Callback(hObject, eventdata, handles)
% hObject    handle to surf_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
surf(handles.current_data);

% --- Executes on button press in mesh_button.
function mesh_button_Callback(hObject, eventdata, handles)
% hObject    handle to mesh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mesh(handles.current_data);

% --- Executes on button press in contour_button.
function contour_button_Callback(hObject, eventdata, handles)
% hObject    handle to contour_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contour(handles.current_data);

% --- Executes on selection change in peaks_popup.
function peaks_popup_Callback(hObject, eventdata, handles)
% hObject    handle to peaks_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: contents = cellstr(get(hObject,'String')) returns peaks_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from peaks_popup
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
case 'peaks' % User selects peaks.
   handles.current_data = handles.peaks;
case 'membrane' % User selects membrane.
   handles.current_data = handles.membrane;
case 'sinc' % User selects sinc.
   handles.current_data = handles.sinc;
end
% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function peaks_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peaks_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
