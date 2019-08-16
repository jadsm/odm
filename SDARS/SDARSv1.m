function varargout = SDARSv1(varargin)
% SDARSV1 MATLAB code for SDARSv1.fig
%      SDARSV1, by itself, creates a new SDARSV1 or raises the existing
%      singleton*.
%
%      H = SDARSV1 returns the handle to a new SDARSV1 or the handle to
%      the existing singleton*.
%
%      SDARSV1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SDARSV1.M with the given input arguments.
%
%      SDARSV1('Property','Value',...) creates a new SDARSV1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SDARSv1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SDARSv1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Copyright
% This code is protected by AstraZeneca's copyright
% The code, however, can be freely distributed, used or modified at will,
% as long as the original publication is correctly cited. 
% The citation should say (or similar):
% "Delgado San Martin et al. (2015)
% Tumour stromal morphology impacts nanomedicine cytotoxicity
% in patient-derived xenografts. Nanomedicine: NBM."
%
% Edit the above text to modify the response to help SDARSv1

% Last Modified by GUIDE v2.5 28-Nov-2014 15:51:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SDARSv1_OpeningFcn, ...
                   'gui_OutputFcn',  @SDARSv1_OutputFcn, ...
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


% --- Executes just before SDARSv1 is made visible.
function SDARSv1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SDARSv1 (see VARARGIN)

% Choose default command line output for SDARSv1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SDARSv1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SDARSv1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 [Macro] = handles2Macro(handles);
 SDARSMain(Macro)


function CD11_Callback(hObject, eventdata, handles)
% hObject    handle to CD11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD11 as text
%        str2double(get(hObject,'String')) returns contents of CD11 as a double


% --- Executes during object creation, after setting all properties.
function CD11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD3_Callback(hObject, eventdata, handles)
% hObject    handle to CD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD3 as text
%        str2double(get(hObject,'String')) returns contents of CD3 as a double


% --- Executes during object creation, after setting all properties.
function CD3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD2_Callback(hObject, eventdata, handles)
% hObject    handle to CD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD2 as text
%        str2double(get(hObject,'String')) returns contents of CD2 as a double


% --- Executes during object creation, after setting all properties.
function CD2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD1_Callback(hObject, eventdata, handles)
% hObject    handle to CD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD1 as text
%        str2double(get(hObject,'String')) returns contents of CD1 as a double


% --- Executes during object creation, after setting all properties.
function CD1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD6_Callback(hObject, eventdata, handles)
% hObject    handle to CD6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD6 as text
%        str2double(get(hObject,'String')) returns contents of CD6 as a double


% --- Executes during object creation, after setting all properties.
function CD6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD5_Callback(hObject, eventdata, handles)
% hObject    handle to CD5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD5 as text
%        str2double(get(hObject,'String')) returns contents of CD5 as a double


% --- Executes during object creation, after setting all properties.
function CD5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD4_Callback(hObject, eventdata, handles)
% hObject    handle to CD4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD4 as text
%        str2double(get(hObject,'String')) returns contents of CD4 as a double


% --- Executes during object creation, after setting all properties.
function CD4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD9_Callback(hObject, eventdata, handles)
% hObject    handle to CD9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD9 as text
%        str2double(get(hObject,'String')) returns contents of CD9 as a double


% --- Executes during object creation, after setting all properties.
function CD9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD8_Callback(hObject, eventdata, handles)
% hObject    handle to CD8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD8 as text
%        str2double(get(hObject,'String')) returns contents of CD8 as a double


% --- Executes during object creation, after setting all properties.
function CD8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CD7_Callback(hObject, eventdata, handles)
% hObject    handle to CD7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD7 as text
%        str2double(get(hObject,'String')) returns contents of CD7 as a double


% --- Executes during object creation, after setting all properties.
function CD7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CD12.
function CD12_Callback(hObject, eventdata, handles)
% hObject    handle to CD12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CD12



function CD10_Callback(hObject, eventdata, handles)
% hObject    handle to CD10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CD10 as text
%        str2double(get(hObject,'String')) returns contents of CD10 as a double


% --- Executes during object creation, after setting all properties.
function CD10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CD10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AR1_Callback(hObject, eventdata, handles)
% hObject    handle to AR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AR1 as text
%        str2double(get(hObject,'String')) returns contents of AR1 as a double


% --- Executes during object creation, after setting all properties.
function AR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AR2_Callback(hObject, eventdata, handles)
% hObject    handle to AR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AR2 as text
%        str2double(get(hObject,'String')) returns contents of AR2 as a double


% --- Executes during object creation, after setting all properties.
function AR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N1_Callback(hObject, eventdata, handles)
% hObject    handle to N1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N1 as text
%        str2double(get(hObject,'String')) returns contents of N1 as a double


% --- Executes during object creation, after setting all properties.
function N1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EO1.
function EO1_Callback(hObject, eventdata, handles)
% hObject    handle to EO1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EO1


% --- Executes on button press in EO2.
function EO2_Callback(hObject, eventdata, handles)
% hObject    handle to EO2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EO2


% --- Executes on button press in EO3.
function EO3_Callback(hObject, eventdata, handles)
% hObject    handle to EO3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EO3



function ST3_Callback(hObject, eventdata, handles)
% hObject    handle to ST3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ST3 as text
%        str2double(get(hObject,'String')) returns contents of ST3 as a double


% --- Executes during object creation, after setting all properties.
function ST3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ST3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ST2_Callback(hObject, eventdata, handles)
% hObject    handle to ST2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ST2 as text
%        str2double(get(hObject,'String')) returns contents of ST2 as a double


% --- Executes during object creation, after setting all properties.
function ST2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ST2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ST1_Callback(hObject, eventdata, handles)
% hObject    handle to ST1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ST1 as text
%        str2double(get(hObject,'String')) returns contents of ST1 as a double


% --- Executes during object creation, after setting all properties.
function ST1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ST1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ST6_Callback(hObject, eventdata, handles)
% hObject    handle to ST6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ST6 as text
%        str2double(get(hObject,'String')) returns contents of ST6 as a double


% --- Executes during object creation, after setting all properties.
function ST6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ST6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ST5_Callback(hObject, eventdata, handles)
% hObject    handle to ST5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ST5 as text
%        str2double(get(hObject,'String')) returns contents of ST5 as a double


% --- Executes during object creation, after setting all properties.
function ST5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ST5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ST4_Callback(hObject, eventdata, handles)
% hObject    handle to ST4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ST4 as text
%        str2double(get(hObject,'String')) returns contents of ST4 as a double


% --- Executes during object creation, after setting all properties.
function ST4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ST4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PP1_Callback(hObject, eventdata, handles)
% hObject    handle to PP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PP1 as text
%        str2double(get(hObject,'String')) returns contents of PP1 as a double


% --- Executes during object creation, after setting all properties.
function PP1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PP1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PP2_Callback(hObject, eventdata, handles)
% hObject    handle to PP2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PP2 as text
%        str2double(get(hObject,'String')) returns contents of PP2 as a double


% --- Executes during object creation, after setting all properties.
function PP2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PP2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PP3_Callback(hObject, eventdata, handles)
% hObject    handle to PP3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PP3 as text
%        str2double(get(hObject,'String')) returns contents of PP3 as a double


% --- Executes during object creation, after setting all properties.
function PP3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PP3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveMacro.
function SaveMacro_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMacro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadMacro.
function LoadMacro_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMacro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[MacroName,~] = uigetfile([cd '\Macros']);
% Load macro
% decodify handles
load(MacroName);
SDARSv1_OpeningFcn(hObject, eventdata, handles)

% --- Executes on button press in Searchbutton.
function Searchbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Searchbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.FileName handles.PathName] = uigetfile;
% Load image
I = imread([handles.PathName handles.FileName]);
handles.I = I;
% Plot image
gca = handles.axes1;imshow(I)
% re-write label
dotpos = strfind(handles.FileName,'.');
set(handles.LoadImageName,'String',handles.FileName(1:dotpos(end)-1));

SDARSv1_OpeningFcn(hObject, eventdata, handles)
function [handles] = Macro2handles(Macro,handles)

function [Macro] = handles2Macro(handles)
aa =[77,75,65,64,63,61,38,22,23,24,19,20,21,...
    16,17,18,12,25,36,35,56,57,...
    58,53,54,55,49,48,47,46,45,44,41,40,39];
Macro = struct2cell(handles);
Macro = Macro(aa);
for i = 3:length(Macro)
    if i > 26
    Field = 'Value';
    Macro(i) = {get(Macro{i},Field)};
    elseif i==20||i == 24||i == 25||i == 26
    Field = 'String'; 
    Macro(i) = {get(Macro{i},Field)};
    else
    Field = 'String'; 
    Macro(i) = {str2num(get(Macro{i},Field))};
    end
end



function PP4_Callback(hObject, eventdata, handles)
% hObject    handle to PP4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PP4 as text
%        str2double(get(hObject,'String')) returns contents of PP4 as a double


% --- Executes during object creation, after setting all properties.
function PP4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PP4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ST7.
function ST7_Callback(hObject, eventdata, handles)
% hObject    handle to ST7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ST7


% --- Executes on button press in ST8.
function ST8_Callback(hObject, eventdata, handles)
% hObject    handle to ST8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ST8


% --- Executes on button press in ST9.
function ST9_Callback(hObject, eventdata, handles)
% hObject    handle to ST9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ST9


% --- Executes on button press in ST10.
function ST10_Callback(hObject, eventdata, handles)
% hObject    handle to ST10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ST10


% --- Executes on button press in ST11.
function ST11_Callback(hObject, eventdata, handles)
% hObject    handle to ST11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ST11


% --- Executes on button press in ST12.
function ST12_Callback(hObject, eventdata, handles)
% hObject    handle to ST12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ST12
