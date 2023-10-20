function varargout = TrackReanalyzerGui(varargin)
% TRACKREANALYZERGUI M-file for TrackReanalyzerGui.fig
%      TRACKREANALYZERGUI, by itself, creates a new TRACKREANALYZERGUI or raises the existing
%      singleton*.
%
%      H = TRACKREANALYZERGUI returns the handle to a new TRACKREANALYZERGUI or the handle to
%      the existing singleton*.
%
%      TRACKREANALYZERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKREANALYZERGUI.M with the given input arguments.
%
%      TRACKREANALYZERGUI('Property','Value',...) creates a new TRACKREANALYZERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrackReanalyzerGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrackReanalyzerGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrackReanalyzerGui

% Last Modified by GUIDE v2.5 09-Mar-2010 12:47:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrackReanalyzerGui_OpeningFcn, ...
                   'gui_OutputFcn',  @TrackReanalyzerGui_OutputFcn, ...
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


% --- Executes just before TrackReanalyzerGui is made visible.
function TrackReanalyzerGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrackReanalyzerGui (see VARARGIN)

% Choose default command line output for TrackReanalyzerGui
%handles.output = hObject;

handles.mra = MaggotReAnalyzer;
handles.hideAnalyze = false;
handles = assignToStruct (handles, varargin);
handles.currentInd = handles.mra.track.pt(1).ind;
if (handles.hideAnalyze)
    set (handles.rean, 'Visible', 'off','Enable','off');
end

% Update handles structure
guidata(hObject, handles);
updateAll(hObject, handles);
% UIWAIT makes TrackReanalyzerGui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrackReanalyzerGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if (~isempty(handles))
    varargout{1} = handles.mra;
    close(handles.figure1);
else
    varargout{1} = [];
end


% --- Executes on slider movement.
function targetAreaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to targetAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
if (v ~= handles.mra.targetArea)
    handles.mra.targetArea = v;
    guidata(hObject, handles);
    updateAll(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function targetAreaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', 10);
set(hObject, 'Max', 510);
set(hObject, 'SliderStep', [0.002 0.2]);


% --- Executes on slider movement.
function frameNbrSlider_Callback(hObject, eventdata, handles)
% hObject    handle to frameNbrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = round(get(hObject,'Value'));
if (v ~= handles.currentInd)
    handles.currentInd = v;
    guidata(hObject, handles);
    updateAll(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function frameNbrSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameNbrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function frameNbr_Callback(hObject, eventdata, handles)
% hObject    handle to frameNbr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameNbr as text
%        str2double(get(hObject,'String')) returns contents of frameNbr as a double
v = round(str2double(get(hObject,'String')));
if (v ~= handles.currentInd)
    handles.currentInd = v;
    guidata(hObject, handles);
    updateAll(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function frameNbr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameNbr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function targetArea_Callback(hObject, eventdata, handles)
% hObject    handle to targetArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetArea as text
%        str2double(get(hObject,'String')) returns contents of targetArea as a double
v = str2double(get(hObject,'String'));
if (v ~= handles.mra.targetArea)
    handles.mra.targetArea = v;
    guidata(hObject, handles);
    updateAll(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function targetArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function contour_scale_Callback(hObject, eventdata, handles)
% hObject    handle to contour_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contour_scale as text
%        str2double(get(hObject,'String')) returns contents of contour_scale as a double
v = str2double(get(hObject,'String'));
if (v ~= handles.mra.contourScale)
    handles.mra.contourScale = v;
    guidata(hObject, handles);
    updateAll(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function contour_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contour_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_ctr_angle_Callback(hObject, eventdata, handles)
% hObject    handle to max_ctr_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_ctr_angle as text
%        str2double(get(hObject,'String')) returns contents of max_ctr_angle as a double
v = str2double(get(hObject,'String'));
if (v ~= handles.mra.maxContourAngle)
    handles.mra.maxContourAngle = v;
    guidata(hObject, handles);
    updateAll(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function max_ctr_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_ctr_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in retvals.
function retvals_Callback(hObject, eventdata, handles)
% hObject    handle to retvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

% --- Executes on button press in rean.
function rean_Callback(hObject, eventdata, handles)
% hObject    handle to rean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mra.reExtractTrack();
uiresume(handles.figure1);

function updateAll(hObject, handles)
%function updateAll(hObject, handles)
try
    if (~isempty(handles.mra.track) && ~isempty(handles.mra.track.pt))
        mn = handles.mra.track.pt(1).ind;
        mx = handles.mra.track.pt(end).ind;
        set(handles.frameNbrSlider, 'Min', mn);
        set(handles.frameNbrSlider, 'Max', mx);
        if (handles.currentInd < mn)
            handles.currentInd = mn;
            guidata(hObject, handles)
        end
        if (handles.currentInd > mx)
            handles.currentInd = mx;
            guidata(hObject, handles)
        end
        set (handles.frameNbrSlider, 'SliderStep', [1 10] / double(mx-mn));
    end
    set(handles.frameNbrSlider, 'Value', handles.currentInd);
    set(handles.frameNbr, 'String', num2str(handles.currentInd));
    set(handles.contour_scale, 'String', num2str(handles.mra.contourScale));
    set(handles.max_ctr_angle, 'String', num2str(handles.mra.maxContourAngle));
    set(handles.targetAreaSlider, 'Value', handles.mra.targetArea);
    set(handles.targetArea, 'String', num2str(handles.mra.targetArea));
    if (~isempty(handles.mra.track) && ~isempty(handles.mra.track.pt))
        ind = find([handles.mra.track.pt.ind] == handles.currentInd,1,'first');
        if (~isempty(ind))
            pt2 = handles.mra.rethreshold(handles.mra.track.pt(ind));
            pt2 = handles.mra.findHT(pt2);
            axes(handles.axes1);
            pt2.drawTrackImage();
            set(handles.axes1, 'XTick', [], 'YTick', []);
        end
    end
catch me
    disp('ooops.  you did something bad!');
end


