function varargout = bi_chromatic_GUI(varargin)
% BI_CHROMATIC_GUI MATLAB code for bi_chromatic_GUI.fig
%      BI_CHROMATIC_GUI, by itself, creates a new BI_CHROMATIC_GUI or raises the existing
%      singleton*.
%
%      H = BI_CHROMATIC_GUI returns the handle to a new BI_CHROMATIC_GUI or the handle to
%      the existing singleton*.
%
%      BI_CHROMATIC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BI_CHROMATIC_GUI.M with the given input arguments.
%
%      BI_CHROMATIC_GUI('Property','Value',...) creates a new BI_CHROMATIC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bi_chromatic_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bi_chromatic_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bi_chromatic_GUI

% Last Modified by GUIDE v2.5 02-Dec-2014 11:12:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bi_chromatic_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @bi_chromatic_GUI_OutputFcn, ...
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


% --- Executes just before bi_chromatic_GUI is made visible.
function bi_chromatic_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bi_chromatic_GUI (see VARARGIN)

% Choose default command line output for bi_chromatic_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bi_chromatic_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bi_chromatic_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in StartAnalysis.
function StartAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to StartAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k = str2num(get(handles.k, 'String'));
Delta = str2num(get(handles.Delta, 'String'));
d_1 = str2num(get(handles.d_1Val, 'String'));
d_2 = str2num(get(handles.d_2Val, 'String'));
Z = str2num(get(handles.Z, 'String'));
V = str2num(get(handles.V, 'String'));
[Vals] = floquet_matrix(k, Delta, d_1, d_2, Z, V);


A = Vals(1,1);
B = Vals(2,1);
C = Vals(3,1);
D = Vals(4,1);
E = Vals(5,1);

t = 100e-6;
npts = 1e6;
rho_12 = zeros(1, npts);

for i = (0:npts),
    time = i*t;
    phi_1 = d_1*time;
    phi_2 = (d_1 + d_2)*time;
    rho = A+B*exp(1i*phi_1)+C*exp(1i*phi_2)+D*exp(-1i*phi_1)+E*exp(-1i*phi_2);
    rho_12(i+1) = real(rho);
end;
plot(rho_12(1:npts));

function k_Callback(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k as text
%        str2double(get(hObject,'String')) returns contents of k as a double


% --- Executes during object creation, after setting all properties.
function k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Delta_Callback(hObject, eventdata, handles)
% hObject    handle to Delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Delta as text
%        str2double(get(hObject,'String')) returns contents of Delta as a double


% --- Executes during object creation, after setting all properties.
function Delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_1Val_Callback(hObject, eventdata, handles)
% hObject    handle to d_1Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_1Val as text
%        str2double(get(hObject,'String')) returns contents of d_1Val as a double


% --- Executes during object creation, after setting all properties.
function d_1Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_1Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_2Val_Callback(hObject, eventdata, handles)
% hObject    handle to d_2Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_2Val as text
%        str2double(get(hObject,'String')) returns contents of d_2Val as a double


% --- Executes during object creation, after setting all properties.
function d_2Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_2Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Z_Callback(hObject, eventdata, handles)
% hObject    handle to Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z as text
%        str2double(get(hObject,'String')) returns contents of Z as a double


% --- Executes during object creation, after setting all properties.
function Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function V_Callback(hObject, eventdata, handles)
% hObject    handle to V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of V as text
%        str2double(get(hObject,'String')) returns contents of V as a double


% --- Executes during object creation, after setting all properties.
function V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
