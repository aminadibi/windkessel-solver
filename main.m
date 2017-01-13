%************************* Windkessel Model Solver ***********************
%*******************       Coded by Amin Adibi     ***********************
%       As a Project for Fund. of Electrical Engineering Course 
% 
%
%  Student No. 86100332
%  Contact: ma.adibi@gmail.com
%
%  January 2012 - Tehran
%  Sharif University of Technology
%
%*************************************************************************
%*************************************************************************




function varargout = Windkessel(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 02-Jan-2012 12:20:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Loading circuit images into GUI 

backgroundImage = importdata('2D.png');
%select the axes
axes(handles.axes_2D);
%place image onto the axes
image(backgroundImage);
%remove the axis tick marks
axis off

backgroundImage = importdata('3D.png');
%select the axes
axes(handles.axes_3D);
%place image onto the axes
image(backgroundImage);
%remove the axis tick marks
axis off


backgroundImage = importdata('4D.png');
%select the axes
axes(handles.axes_4D);
%place image onto the axes
image(backgroundImage);
%remove the axis tick marks
axis off
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.

function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radio_2el.
function radio_2el_Callback(hObject, eventdata, handles)
% hObject    handle to radio_2el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_2el


% --------------------------------------------------------------------
function Menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radio_3el.
function radio_3el_Callback(hObject, eventdata, handles)
% hObject    handle to radio_3el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_3el


% --- Executes on button press in radio_4el.
function radio_4el_Callback(hObject, eventdata, handles)
% hObject    handle to radio_4el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_4el

function Parameter_R_Callback(hObject, eventdata, handles)
% hObject    handle to Parameter_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Parameter_R as text
%        str2double(get(hObject,'String')) returns contents of Parameter_R as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Parameter_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Parameter_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Parameter_C_Callback(hObject, eventdata, handles)
% hObject    handle to Parameter_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Parameter_C as text
%        str2double(get(hObject,'String')) returns contents of Parameter_C as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Parameter_C_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Parameter_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Parameter_r_Callback(hObject, eventdata, handles)
% hObject    handle to Parameter_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Parameter_r as text
%        str2double(get(hObject,'String')) returns contents of Parameter_r as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Parameter_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Parameter_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Parameter_L_Callback(hObject, eventdata, handles)
% hObject    handle to Parameter_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Parameter_L as text
%        str2double(get(hObject,'String')) returns contents of Parameter_L as a double
input = str2num(get(hObject,'String'));

%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Parameter_L_CreateFcn(hObject, ~, handles)
% hObject    handle to Parameter_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
function res1 = eq3(t,I,r,R,C,h) % ODE number 2
    res1=(-(r+R)/(r*R*C)).*I+ p_mom(t,h)./(r*R*C) + dp_mom(t,h)./r;

% The Following two ODEs describe 4-Element Windkessel Model
function res1 = eq41(I2,I1,index,h) % ODE number 4 - part 1       
        if ((index-1)<1e-4) res1=0.01; %  Handles a special case in the numerical calculation of dp/dt
        else      
            res1 = (I2-I1)./h; 
        end
            
function z_prime = eq42(t,I,z,r,R,C,L,h) % ODE number 4 - part 2, z means dIdt
z_prime=(-(r*R*C+C)/(R*L*C)).*z - ((R+r)./(L*C*R)).*I + p_mom(t,h)./(L*R*C) + dp_mom(t,h)./L;    
    
function p_res = p_mom(t,h) % Pressure
    if t==-h p_res=83.6017 % Handles a special case in the numerical calculation of dp/dt
    else
        p_res = (13.128-0.423*cos(2*pi*t)+1.417*sin(2*pi*t)-0.418*cos(4*pi*t)+0.747*sin(4*pi*t)-0.615*cos(6*pi*t)+0.312*sin(6*pi*t)-0.326*cos(8*pi*t)-0.228*sin(8*pi*t)+0.068*cos(10*pi*t)-0.0506*sin(10*pi*t)-0.082*cos(14*pi*t)+0.093*sin(14*pi*t)-0.111*cos(16*pi*t)-0.0834*sin(16*pi*t)+0.035*cos(18*pi*t)-0.015*sin(18*pi*t)-0.074*cos(20*pi*t)+0.028*sin(20*pi*t)-0.036*cos(22*pi*t)-0.048*sin(22*pi*t))*7.5006;
    end
 

function dp_res = dp_mom(t,h) % Numerical function for dp/dt
    dp_res=(p_mom(t,h)-p_mom(t-h,h))./h;
    
    
% --- Executes on button press in pushbutton_simulate.
function pushbutton_simulate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Here we read parameters from GUI
R_text = get(handles.Parameter_R,'String');
R = str2num (R_text);

C_text = get(handles.Parameter_C,'String');
C = str2num (C_text);

r_text = get(handles.Parameter_r,'String');
r = str2num (r_text);

L_text = get(handles.Parameter_L,'String');
L = str2num (L_text);

T_text = get(handles.input_T,'String');
T = str2num (T_text);

h_text = get(handles.input_h,'String');
h = str2num (h_text);
% End of "Reading Parameters"

I_temp(1) = 0.01; %Initial Value

% Beginning of Solving Windkessel Model
index = 1;     % Used to keep track in RK4 for loop
I(1)= 0.01;    % Blood flow rate, initial condition

I_temp(1) = 0.01;

 if get(handles.radio_2el,'Value') == 1 % The case for 2 Elements model
     % Here we have an analytical equation that gives flow rate (I)
     % based on known Pressure Wave and its derviative
     for t = 0 : h : T 
         if ((T-t) < 1e-4), break, end % handles a special case when h<<1
         I (index+1) = p_mom(t,h)./(R) + C .* dp_mom(t,h);
         index = index + 1;
     end
 end
 
     
 if get(handles.radio_3el,'Value') == 1 % The case for 3 Elements model
     % Here we've got an ODE to solve. We use classic 4-parameter
     % Runge-Kutta method to solve it.
     
    for t = 0 : h : T 
        if ((T-t) < 1e-4), break, end % handles a special case when h<<1
    
        k11 = eq3(t,I_temp(1),r,R,C,h);
        I_temp(2) = I_temp(1) + k11*h/2;
    
        k21 = eq3(t+h/2,I_temp(2),r,R,C,h);
        I_temp(2) = I_temp(1) + k21*h/2;    
    
        k31 = eq3(t+h/2,I_temp(2),r,R,C,h);
        I_temp(3) = I_temp(1) + k31*h;
    
        k41= eq3(t+h,I_temp(3),r,R,C,h);
    
        k1 = 1/6 * (k11 + 2*k21 + 2*k31 + k41);
    
        I_temp(3) = I_temp(1) + k1*h;
    
        I(index+1) = I(index) + k1*h;  %Assigning new values
    
        index = index + 1;
    
        I_temp = [I_temp(3), 0, 0];   % reseting temporary values
    
    end
 end

 if get(handles.radio_4el,'Value') == 1 % The case for 4 Elements model
     
    % Here we have a second order ODE.
    % We will convert it into a two first order ODEs and solve them 
    % simultanously using Classic 4-Parameter Runge-Kutta Method
     I (2) = 0.011; %Extra Initial Values 
     I_temp(2) = 0.011;
     index=1;
     dIdt=0;
     for t = 0 : h : T 
         if ((T-t) < 1e-4), break, end % handles a special case when h<<1
         
        k11 = eq41(I_temp(1), I(index),index,h);
        dIdt(2) = dIdt(1) + k11*h/2;
        
        k21 = eq42(t,I_temp(1),dIdt(2),r,R,C,L,h);
        I_temp(2) = I_temp(1) + k21*h/2;
       
        k12 = eq41 (I_temp(2), I(index),index,h);
        dIdt (2) = dIdt(1) + k12*h/2;
    
        k22 = eq42(t+h/2,I_temp(2),dIdt(2),r,R,C,L,h);
        I_temp(2) = I_temp(1) + k22*h/2;
        
        k13 = eq41(I_temp(2), I(index),index,h);
        dIdt(3) = dIdt(1) + k13*h;
    
        k23 = eq42(t+h/2,I_temp(2),dIdt(3),r,R,C,L,h);
        I_temp(3) = I_temp(1) + k23*h;
        
        k14= eq41(I_temp(3), I(index),index,h);
        k1 = 1/6 * (k11 + 2*k12 + 2*k13 + k14);
        
        k24= eq42(t+h,I_temp(3),dIdt(3),r,R,C,L,h);
        k2 = 1/6 * (k21 + 2*k22 + 2*k23 + k24);

        dIdt(3) = dIdt(1) + k1*h;
        I_temp(3) = I_temp(1) + k2*h;
    
        I(index+1) = I(index) + k2*h;  %Assigning new values
    
        index = index + 1;
        dIdt = [dIdt(3),0,0];
        I_temp = [I_temp(3), 0, 0];   % reseting temporary values

      end
 end 

% Start Drawing Plot
Time = 0 : h : T;

%plots the x and y data
plot(handles.Output_Graph,Time,I);

%Feshar = (13.128-0.423*cos(2*pi*Time)+1.417*sin(2*pi*Time)-0.418*cos(4*pi*Time)+0.747*sin(4*pi*Time)-0.615*cos(6*pi*Time)+0.312*sin(6*pi*Time)-0.326*cos(8*pi*Time)-0.228*sin(8*pi*Time)+0.068*cos(10*pi*Time)-0.0506*sin(10*pi*Time)-0.082*cos(14*pi*Time)+0.093*sin(14*pi*Time)-0.111*cos(16*pi*Time)-0.0834*sin(16*pi*Time)+0.035*cos(18*pi*Time)-0.015*sin(18*pi*Time)-0.074*cos(20*pi*Time)+0.028*sin(20*pi*Time)-0.036*cos(22*pi*Time)-0.048*sin(22*pi*Time))*7.5006;

%adds a title, x-axis description, and y-axis description
%title('Volumetric Flow Rate');
xlabel('Time');
ylabel('Flow Rate');
%guidata(hObject, handles); %updates the handles

% End of Solving Windkessel Model
guidata(hObject, handles);


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%clear handles.Output_Graph;
%clf (handles.Output_Graph);
guidata(hObject, handles);


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


% --- Executes on key press with focus on pushbutton_simulate and none of its controls.
function pushbutton_simulate_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_simulate (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function input_T_Callback(hObject, eventdata, handles)
% hObject    handle to input_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_T as text
%        str2double(get(hObject,'String')) returns contents of input_T as a double


% --- Executes during object creation, after setting all properties.
function input_T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function input_h_Callback(hObject, eventdata, handles)
% hObject    handle to input_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_h as text
%        str2double(get(hObject,'String')) returns contents of input_h as a double


% --- Executes during object creation, after setting all properties.
function input_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_h (see GCBO)
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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
