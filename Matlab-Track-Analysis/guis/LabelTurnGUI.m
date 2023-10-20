function varargout = LabelTurnGUI(varargin)
% LABELTURNGUI M-file for LabelTurnGUI.fig
%      LABELTURNGUI, by itself, creates a new LABELTURNGUI or raises the existing
%      singleton*.
%
%      H = LABELTURNGUI returns the handle to a new LABELTURNGUI or the handle to
%      the existing singleton*.
%
%      LABELTURNGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LABELTURNGUI.M with the given input arguments.
%
%      LABELTURNGUI('Property','Value',...) creates a new LABELTURNGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LabelTurnGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LabelTurnGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LabelTurnGUI

% Last Modified by GUIDE v2.5 17-Sep-2010 15:22:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LabelTurnGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LabelTurnGUI_OutputFcn, ...
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


% --- Executes just before LabelTurnGUI is made visible.
function LabelTurnGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LabelTurnGUI (see VARARGIN)

% Choose default command line output for LabelTurnGUI
handles.output = hObject;

if (isempty(varargin))
    error ('no args passed. TrackPlayer(src, options), where src is eset(s), expt(s), track(s)');
end
src = varargin{1};
%tic
global ltg_alltracks;
if (isa (src, 'ExperimentSet'))
    e = [src.expt];
    ltg_alltracks= [e.track];
else
    if (isa (src, 'Experiment'))
        ltg_alltracks = [src.track];
    else
        if (isa (src, 'Track'))
            ltg_alltracks = src;
        else
            error ('LabelTurnGui(src, options), where src is eset(s), expt(s), track(s)');
        end
    end
end
%toc
handles.track = ltg_alltracks(1);
handles.type1 = 'speed';
handles.type2 = 'path';
set(handles.tnumBox, 'String', cellfun(@num2str, num2cell(1:length(ltg_alltracks)),'UniformOutput', false));
topfields = {'path', 'speed', 'deltatheta', 'covRatio'};
fieldList = [topfields, setdiff(unique([fieldnames(handles.track.dq); (handles.track.validDQName)']), topfields)];
set(handles.data1Type, 'String', fieldList);
set(handles.data2Type, 'String', fieldList);
set(handles.data1Type, 'Value', 2);
set(handles.data1Type, 'Value', 1);
handles.delayms = 100;
handles.ptbuffer = 6;
handles = newTrack(hObject, handles, 1);
%toc
% Update handles structure
guidata(hObject, handles);
%toc
% UIWAIT makes LabelTurnGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LabelTurnGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in omegaButton.
function omegaButton_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% hObject    handle to omegaButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of omegaButton
handles = setTurnType (hObject, handles, 'O');
guidata(hObject, handles);

% --- Executes on button press in reverseButton.
function reverseButton_Callback(hObject, eventdata, handles)
% hObject    handle to reverseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reverseButton
handles = setTurnType (hObject, handles, 'R');
guidata(hObject, handles);

% --- Executes on button press in doubleButton.
function doubleButton_Callback(hObject, eventdata, handles)
% hObject    handle to doubleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of doubleButton
handles = setTurnType (hObject, handles, 'D');
guidata(hObject, handles);

% --- Executes on button press in multiButton.
function multiButton_Callback(hObject, eventdata, handles)
% hObject    handle to multiButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiButton
handles = setTurnType (hObject, handles, 'M');
guidata(hObject, handles);

% --- Executes on button press in pauseButton.
function pauseButton_Callback(hObject, eventdata, handles)
% hObject    handle to pauseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pauseButton
handles = setTurnType (hObject, handles, 'P');
guidata(hObject, handles);

% --- Executes on button press in noturnButton.
function noturnButton_Callback(hObject, eventdata, handles)
% hObject    handle to noturnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noturnButton
handles = setTurnType (hObject, handles, 'N');
guidata(hObject, handles);

% --- Executes on button press in noideaButton.
function noideaButton_Callback(hObject, eventdata, handles)
% hObject    handle to noideaButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noideaButton
handles = setTurnType (hObject, handles, 'C');
guidata(hObject, handles);

% --- Executes on selection change in data1Type.
function data1Type_Callback(hObject, eventdata, handles)
% hObject    handle to data1Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data1Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data1Type
contents = cellstr(get(hObject,'String')) ;
handles.type1 =  contents{get(hObject,'Value')} ;
handles = plotDataFields(hObject, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function data1Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data1Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in data2Type.
function data2Type_Callback(hObject, eventdata, handles)
% hObject    handle to data2Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data2Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        data2Type

contents = cellstr(get(hObject,'String')) ;
handles.type2 =  contents{get(hObject,'Value')} ;
handles = plotDataFields(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function data2Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data2Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in playMovie.
function playMovie_Callback(hObject, eventdata, handles)
% hObject    handle to playMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = playMovie(hObject, handles);
guidata(hObject, handles);

function delaymsedit_Callback(hObject, eventdata, handles)
% hObject    handle to delaymsedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delaymsedit as text
%        str2double(get(hObject,'String')) returns contents of delaymsedit as a double
handles.delayms = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function delaymsedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delaymsedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ptbufferEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ptbufferEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptbufferEdit as text
%        str2double(get(hObject,'String')) returns contents of ptbufferEdit as a double
handles.ptbuffer =  str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ptbufferEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptbufferEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prevButton.
function prevButton_Callback(hObject, eventdata, handles)
% hObject    handle to prevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = newSharpTurn(hObject, handles, handles.sti - 1);
guidata(hObject, handles);

% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = newSharpTurn(hObject, handles, handles.sti + 1);
guidata(hObject, handles);


% --- Executes on selection change in tnumBox.
function tnumBox_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to tnumBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tnumBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tnumBox
handles = newTrack(hObject, handles, get(hObject, 'Value'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tnumBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tnumBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = newTrack(hObject, handles, tracknum)
global ltg_alltracks;
tracknum = max(1, tracknum);
tracknum = min(tracknum, length(ltg_alltracks));

if (isfield(handles, 'tnum') && handles.tnum == tracknum)
    return;
end

handles.tnum = tracknum;
handles.track = ltg_alltracks(handles.tnum);
handles.st = handles.track.sharpTurn;
handles = newSharpTurn(hObject, handles, 1);

function handles = newSharpTurn(hObject, handles, stnum)
stnum = max(1, stnum);
stnum = min(stnum, length(handles.st));
handles.sti = stnum;
if (stnum > 0)
    handles.stc = handles.st(handles.sti);
end
set(handles.compLabel, 'String', handles.stc.type);
set(handles.userLabel, 'String', handles.stc.usertype);
handles = plotDataFields(hObject, handles);
if (stnum > 0 && ~isfinite(handles.stc.userCode))
    playMovie(hObject, handles);
end
guidata(hObject, handles);


function handles = playMovie(hObject, handles)
handles.stc.playMovie('ptbuffer', handles.ptbuffer, 'delayTime', handles.delayms/1000, ...
        'LabelTurns', false, 'Axes', handles.movAxes, 'axisSize', 20);
    
function handles = setTurnType (hObject, handles, ttype)

handles.stc.setUserType(ttype);
handles = newSharpTurn(hObject, handles, handles.sti + 1);

function handles = plotDataFields (hObject, handles)
Axes = [handles.dataAxes1, handles.dataAxes2];
types = {handles.type1, handles.type2};
for j = 1:2
    cla(Axes(j));
    if (strcmpi(types{j},'path'))
        handles.stc.plotPath('sloc', 'k-', 'Axes', Axes(j));
    else 
        handles.stc.plotFields('eti', types{j}, 'Axes', Axes(j));
    end
end
