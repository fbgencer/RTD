function varargout = rtd(varargin)
% RTD MATLAB code for rtd.fig
%Written by Fikret Başar GENCER, FBG.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rtd_OpeningFcn, ...
                   'gui_OutputFcn',  @rtd_OutputFcn, ...
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


% --- Executes just before rtd is made visible.
function rtd_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rtd (see VARARGIN)
datacursormode on
warning off;
clc;
handles.precision = 50;
handles.applied_voltage = 1; %this is iterator for plots
handles.applied_voltage_low = 0;
handles.applied_voltage_high = 0.5;
handles.applied_voltage_size = 10;
handles.wave_energy_low = 0;
handles.wave_energy_high = 0.5;
handles.wave_energy_size = 100;

handles.wave_energy_slider.Min = handles.wave_energy_low;
handles.wave_energy_slider.Max = handles.wave_energy_high;
handles.wave_energy_slider.Value = handles.wave_energy_high;

%transmission,reflection and regions
handles.wave_energy = [];
handles.t = {};
handles.r = {};
handles.region_matrix = {};
handles.current = [];
handles.temp = 300;

handles.simulation_type = "tcoef_ivcurve";

%Table
handles.uitable.Data = [2, 0.5];
handles.uitable.RowName = {'Barrier'};
table_widths = handles.uitable.Data(:,1);
table_potentials = handles.uitable.Data(:,2);
%Plot first axes 
update_structure(hObject, eventdata, handles,table_widths,table_potentials);

%Axes properties
handles.yaxis_type = 'Log';


% Choose default command line output for rtd
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rtd wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = rtd_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
case '1-Barrier'
   handles.uitable.RowName = {'Barrier'};
   handles.uitable.Data = [ transpose([2]), transpose([0.5]) ]; 
case '2-Barrier' % User selects peaks.
   handles.uitable.RowName = {'Barrier'; 'Well'; 'Barrier'};
   handles.uitable.Data = [ transpose([2 5 2]), transpose([0.5 0 0.5]) ];
   %handles.uitable.Data = transpose([2 5 2]);
case '3-Barrier' % User selects peaks.
   handles.uitable.RowName = {'Barrier'; 'Well'; 'Barrier'; 'Well'; 'Barrier'};
   handles.uitable.Data = [ transpose([2 5 2 5 2]), transpose([0.5 0 0.5 0 0.5]) ];
case '4-Barrier' % User selects peaks.
   handles.uitable.RowName = {'Barrier'; 'Well'; 'Barrier'; 'Well'; 'Barrier';'Well'; 'Barrier'};
  handles.uitable.Data = [ transpose([2 5 2 5 2  5 2]), transpose([0.5 0 0.5 0 0.5 0 0.5]) ];
case '5-Barrier' % User selects peaks.
   handles.uitable.RowName = {'Barrier'; 'Well'; 'Barrier'; 'Well'; 'Barrier';'Well'; 'Barrier';'Well'; 'Barrier'};
   handles.uitable.Data = [ transpose([2 5 2 5 2 5 2 5 2]), transpose([0.5 0 0.5 0 0.5 0 0.5 0 0.5]) ];
case 'Custom'
        
end

table_widths = handles.uitable.Data(:,1);
table_potentials = handles.uitable.Data(:,2);

update_structure(hObject, eventdata, handles,table_widths,table_potentials);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in simulation_type_popup.
function simulation_type_popup_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');

switch str{val}
case 'T.Coef + I-V Curve'
    handles.simulation_type = "tcoef_ivcurve";
case 'Transmission Coefficient' % User selects peaks.
    handles.simulation_type = "tcoef";
case 'I-V Curve' % User selects peaks.
    handles.simulation_type = "ivcurve";
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function simulation_type_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simulation_type_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes when entered data in editable cell(s) in uitable.
function uitable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
table_widths = handles.uitable.Data(:,1);
table_potentials = handles.uitable.Data(:,2);
datacursormode on
update_structure(hObject, eventdata, handles,table_widths,table_potentials);
guidata(hObject,handles);

% --- Executes when selected cell(s) is changed in uitable.
function uitable_CellSelectionCallback(hObject, eventdata, handles)

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function precision_text_Callback(hObject, ~, handles)
handles.precision = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function precision_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to precision_text (see GCBO)
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


% --- Executes on selection change in plot_type_popup.
function plot_type_popup_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
disp(handles.simulation_type);
switch str{val}
case 'Transmission Coefficient'
    if(handles.simulation_type == "tcoef" || handles.simulation_type == "tcoef_ivcurve")
        
        handles.plot2 = plot(handles.axes2,handles.wave_energy,handles.t{handles.applied_voltage});
        xlabel(handles.axes2,'energy(eV)')
        grid on;
        ylabel(handles.axes2,'T(E)')
        handles.axes2.YScale = handles.yaxis_type;
    end
case 'Current - Voltage' % User selects peaks.
   if(handles.simulation_type == "ivcurve" || handles.simulation_type == "tcoef_ivcurve")
        %cla(handles.axes2,'reset');
        AV = linspace(handles.applied_voltage_low,handles.applied_voltage_high,handles.applied_voltage_size);
        Temp = handles.temp;
        plot(handles.axes2,AV,handles.current,'b -','LineWidth',1.5);
        ylabel(handles.axes2,'$(J) (A/m^{2})$','Interpreter','latex','FontSize',14,'FontWeight','bold')
        xlabel(handles.axes2,'$AV (V)$','Interpreter','latex','FontSize',14,'FontWeight','bold')
        plot_label = [ 'T ',num2str(Temp),'$^oK$'] ;
        legend(handles.axes2,plot_label);
        handles.axes2.Legend.Location ='northeast';
        handles.axes2.Legend.Interpreter = "Latex";
        handles.axes2.YScale = handles.yaxis_type;
        grid(handles.axes2,'on');
    
   end
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function plot_type_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function yaxis_popup_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
case 'Y-axis: Log'
    handles.yaxis_type = 'Log';
case 'Y-axis: Linear'
    handles.yaxis_type = 'Linear';     
end
handles.axes2.YScale = handles.yaxis_type;
guidata(hObject,handles);


function yaxis_popup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_structure(hObject, eventdata, handles,table_widths,table_potentials)
q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;
hbar =1.0545718e-34; me = 9.110e-31; q_e =1.602e-19;
um = 1e-6; nm = 1e-9;
me = 0.063*me;
%me = 0.0919*me;
kB = 1.38 *1e-23;

left_contact_width = 10;
right_contact_width = 10;

applied_voltage = 0;

potentials = [-applied_voltage table_potentials' 0]*eV; 
widths = [left_contact_width table_widths' right_contact_width]*nm;
potential_profile = @(x) (-applied_voltage*x/( widths(1)-sum(widths(1:end-1)) ) + ...
applied_voltage*(sum(widths(1:end-1))/(widths(1)-sum(widths(1:end-1)) ) ) );

[t,r,region_matrix,k,interface_x] = trans_coef(handles.precision,potentials,widths,0,1,potential_profile);
plot_regions(region_matrix,handles.axes1);
guidata(hObject,handles);


function wave_energy_slider_Callback(hObject, eventdata, handles)
if(get(hObject,'Value') ~= handles.wave_energy_low)
    handles.axes2.XLim = [handles.wave_energy_low,get(hObject,'Value')];
end
guidata(hObject,handles);

function wave_energy_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function wave_energy_low_Callback(hObject, eventdata, handles)
handles.wave_energy_low = str2double(get(hObject,'String'));
handles.wave_energy_slider.Min = handles.wave_energy_low;
handles.wave_energy_slider.Max = handles.wave_energy_high;
guidata(hObject,handles);

function wave_energy_low_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wave_energy_high_Callback(hObject, eventdata, handles)
handles.wave_energy_high = str2double(get(hObject,'String'));
handles.wave_energy_slider.Min = handles.wave_energy_low;
handles.wave_energy_slider.Max = handles.wave_energy_high;
guidata(hObject,handles);

function wave_energy_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simulate_button.
function simulate_button_Callback(hObject, eventdata, handles)
q =1.602e-19;
um = 1e-6; nm = 1e-9;
eV = 1.6*10^-19;
hbar =1.0545718e-34; h = 2*pi*hbar;
me = 9.110e-31;
me = 0.063*me;
%me = 0.0919*me;
kB = 1.38 *1e-23;

left_contact_width = 10;
right_contact_width = 10;

AV = linspace(handles.applied_voltage_low,handles.applied_voltage_high,handles.applied_voltage_size);
widths = [left_contact_width transpose(handles.uitable.Data(:,1)) right_contact_width]*nm;
wave_amplitude = 1;
handles.wave_energy = linspace(handles.wave_energy_low,handles.wave_energy_high,handles.wave_energy_size);

if(handles.simulation_type == "tcoef_ivcurve" || handles.simulation_type == "tcoef")
for iter = 1:size(AV,2)

percent =  round(100*iter/handles.applied_voltage_size);   
handles.simulate_button.String = "Simulating(1/2) " + num2str(percent) + "%";
drawnow;

applied_voltage = AV(1,iter);
potentials = [-applied_voltage transpose(handles.uitable.Data(:,2)) 0]*eV; 

potential_profile = @(x) (-applied_voltage*x/( widths(1)-sum(widths(1:end-1)) ) + ...
applied_voltage*(sum(widths(1:end-1))/(widths(1)-sum(widths(1:end-1)) ) ) );

[handles.t{iter},handles.r{iter},handles.region_matrix{iter},k,interface_x] = ...
trans_coef(handles.precision,potentials,widths,handles.wave_energy*eV,wave_amplitude,potential_profile);
end
end
if(handles.simulation_type == "tcoef_ivcurve" || handles.simulation_type == "ivcurve")
    disp("IV curve running..");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:size(AV,2)
percent =  round(100*iter/handles.applied_voltage_size);   
handles.simulate_button.String = "Simulating(2/2)" + num2str(percent) + "%";
drawnow;
applied_voltage = AV(iter);
potentials = [-applied_voltage transpose(handles.uitable.Data(:,2)) 0]*eV;
potential_profile = @(x) (-applied_voltage*x/( widths(1)-sum(widths(1:end-1)) ) + ...
applied_voltage*(sum(widths(1:end-1))/(widths(1)-sum(widths(1:end-1)) ) ) );
% first we calculate the transmission coefficients of all energies up to
% the fermi energy. transSteps need to be a large number to get all details
% in the transmission spectrum, because the peaks are very sharp
% steps=1;
% EfL=0.005*eV; % fermi energy in left metal
% EfR=0.005*eV; % right metal
% for n=1:steps
%     [t(n),r(n),region_matrix,k,interface_x] = trans_coef(handles.precision,potentials,widths,EfL/steps*(n),wave_amplitude,potential_profile);
%     %T(n)=rtTransmission(AV, EfL/steps*(n), 6);
% end
% % Integrate numerically the expression for the current, integration is
% % performed from 0 to EfL
% I=0;
% for n=1:steps
%     Energy=applied_voltage *q + EfL/steps*(n-1);
%     
%     Fl=DistFermiDirac(Energy,applied_voltage*q+EfL,0);
%     Fr=DistFermiDirac(Energy,EfR,0);
%     I = I + (t(n)*(Fl-Fr)*sqrt(Energy))*EfL/steps;
% end
% handles.current(iter) = (I);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ex burada wave energy cünkü transmission coef'i oyle buluyor
%gerilim altında oluşacak senaryo ise
%J = 4*pi*me*q/(h^3) Integrate[T(Ex)N(Ex)dEx,Emin,Emax];
%numerator = 1+exp(-(Ex-Ef1)/(kB*Temp));
%denominator = 1+exp(-(Ex-Ef2)/(kB*Temp));
%N(Ex) = kB*Temp*log(numerator/denominator);
%AV gerilim altında
%denominator = 1+exp(-(Ex + AV*eV - Ef2)/(kB*Temp));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Emin = 0.005*eV;
Emax = 0.49*eV;
Temp = handles.temp;
Ef = 0.005*eV;
steps = 50;
Ex = linspace(Emin,Emax,steps);

dEx = (Emax-Emin)/steps;
J = 0;
for n=1:steps
    numerator = 1+exp(-(Ex(n)-Ef)/(kB*Temp));
    denominator = 1+exp(-(Ex(n) + applied_voltage*eV - Ef)/(kB*Temp));
    [t,r,region_matrix] = trans_coef(handles.precision,potentials,widths,Ex(n),1,potential_profile);
    N_Ex = log(numerator/denominator);
    %J = 4*pi*me*q/(h^3) Integrate[T(Ex)N(Ex)dEx,Emin,Emax];
    J = J + (t * N_Ex * dEx);
    
%     Energy=applied_voltage * eV + EfL/steps*(n-1);
%     Fl=DistFermiDirac(Energy,applied_voltage * eV+EfL,Temp);
%     Fr=DistFermiDirac(Energy,EfR,Temp);
%     I = I + (t(n)*(Fl-Fr)*sqrt(Energy))*EfL/steps;
end
%I = (q*me*kB*Temp/(2*pi*pi*hbar^3))*I;
J = kB*Temp*4*pi*me*q/(h^3)*J;
handles.current(iter) = (J);
end

end

if(handles.simulation_type == "ivcurve")
    plot(handles.axes2,AV,handles.current,'b -','LineWidth',1.5);
    ylabel(handles.axes2,'$ln(J) (A/m^{2})$','Interpreter','latex','FontSize',14,'FontWeight','bold')
    xlabel(handles.axes2,'$AV (V)$','Interpreter','latex','FontSize',14,'FontWeight','bold')
    plot_label = [ 'T ',num2str(Temp),'$^oK$'] ;
    legend(plot_label);
    ax = gca ;
    ax.Legend.Location ='northeast';
    ax.Legend.Interpreter = "Latex";
    grid on;
    
else 
    plot(handles.axes2,handles.wave_energy,handles.t{handles.applied_voltage},'b -','LineWidth',1.5);
    ylabel(handles.axes2,'$T$','Interpreter','latex','FontSize',14,'FontWeight','bold')
    xlabel(handles.axes2,'$E(eV)$','Interpreter','latex','FontSize',14,'FontWeight','bold')
    grid on;
    %f1 = figure;
    %plot(handles.wave_energy,handles.t{handles.applied_voltage},'b -');
    
    plot_regions(handles.region_matrix{handles.applied_voltage},handles.axes1);
end

handles.axes2.YScale = handles.yaxis_type;
handles.simulate_button.String = 'Simulate'; drawnow;
guidata(hObject,handles);


function applied_voltage_slider_Callback(hObject, eventdata, handles)

if(handles.plot_type_popup.Value == 1 && (handles.simulation_type == "tcoef_ivcurve" || handles.simulation_type == "tcoef") )
handles.applied_voltage = round(map(get(hObject,'Value'),handles.applied_voltage_slider.Min,handles.applied_voltage_slider.Max, ...
    1,handles.applied_voltage_size));

val  =  handles.applied_voltage * (handles.applied_voltage_high - handles.applied_voltage_low)/(handles.applied_voltage_size);
%handles.applied_voltage_low + (get(hObject,'Value')-1)*((handles.applied_voltage_high - handles.applied_voltage_low)/...
%handles.(handles.applied_voltage_slider.Max - handles.applied_voltage_slider.Min));

handles.applied_voltage_text.String = "Applied Voltage:"+num2str(val);


plot(handles.axes2,handles.wave_energy,handles.t{handles.applied_voltage},'b -','LineWidth',1.5);
ylabel(handles.axes2,'$T$','Interpreter','latex','FontSize',14,'FontWeight','bold')
xlabel(handles.axes2,'$E(eV)$','Interpreter','latex','FontSize',14,'FontWeight','bold')
grid on;

plot_regions(handles.region_matrix{handles.applied_voltage},handles.axes1);
handles.axes2.YScale = handles.yaxis_type;
guidata(hObject,handles);
end


function applied_voltage_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function applied_voltage_low_Callback(hObject, eventdata, handles)
handles.applied_voltage_low = str2double(get(hObject,'String'));
guidata(hObject,handles);

function applied_voltage_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function applied_voltage_high_Callback(hObject, eventdata, handles)
handles.applied_voltage_high = str2double(get(hObject,'String'));
guidata(hObject,handles);
function applied_voltage_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function loading_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loading_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function wave_energy_size_text_Callback(hObject, eventdata, handles)
% hObject    handle to wave_energy_size_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wave_energy_size_text as text
handles.wave_energy_size = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function wave_energy_size_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wave_energy_size_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function applied_voltage_size_text_Callback(hObject, eventdata, handles)
handles.applied_voltage_size = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function applied_voltage_size_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to applied_voltage_size_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [res] = map(x, in_min, in_max, out_min,out_max) 
  res = (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;


% --- Executes during object creation, after setting all properties.
function applied_voltage_text_CreateFcn(hObject, eventdata, handles)

function [probability] = DistFermiDirac(E, Ef, T)
    kB = 1.38e-23;
    if(T == 0)
        if(E > Ef)
            probability = 0;
            return;
        else
            probability = 1;
            return;
        end
    else
        denom=exp((E-Ef)/(kB*T)) + 1;
    end
    probability = 1/denom;



function temperature_value_Callback(hObject, eventdata, handles)
handles.temp = str2double(get(hObject,'String')) 


% --- Executes during object creation, after setting all properties.
function temperature_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to temperature_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
