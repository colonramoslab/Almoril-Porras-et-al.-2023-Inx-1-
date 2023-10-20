function varargout = sminspector(varargin)
% SMINSPECTOR M-file for sminspector.fig
%      SMINSPECTOR, by itself, creates a new SMINSPECTOR or raises the existing
%      singleton*.
%
%      H = SMINSPECTOR returns the handle to a new SMINSPECTOR or the handle to
%      the existing singleton*.
%
%      SMINSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMINSPECTOR.M with the given input arguments.
%
%      SMINSPECTOR('Property','Value',...) creates a new SMINSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sminspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sminspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sminspector

% Last Modified by GUIDE v2.5 08-Jun-2010 17:18:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sminspector_OpeningFcn, ...
                   'gui_OutputFcn',  @sminspector_OutputFcn, ...
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


% --- Executes just before sminspector is made visible.
function sminspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sminspector (see VARARGIN)

% Choose default command line output for sminspector
handles.output = hObject;

if (isempty(varargin))
    disp ('sminspector(segmentationModel, pointsource(optional))');
    close(hObject);
    return;
end

global SMI_SM;
global SMI_TRACKS;
SMI_SM = varargin{1};
sc = SMI_SM.segmentationClusters;
if (length(varargin) < 2)    
    kp = [sc.knownPoints];
    SMI_TRACKS = unique([kp.track]);
else
    if (isa(varargin{2}, 'Track'))
        SMI_TRACKS = varargin{2};
    else
        if (isa(varargin{2}, 'Experiment'))
            SMI_TRACKS = [varargin{2}.track];
        else
            if (isa(varargin{2}, 'ExperimentSet'))
                e = [varargin{2}.expt];    
                SMI_TRACKS = [e.track];
            else
                disp ('second argument must be experiment set, experiment, or track');
                close (hObject);
                return;
            end
        end
    end
end
s = cellfun(@num2str, num2cell(1:length(SMI_TRACKS)),'UniformOutput', false);
s = [{'All'} s];
set(handles.tracknumberpopup, 'String',s);
set(handles.xtypepopup, 'String', sc(1).datafields);
set(handles.ytypepopup, 'String', sc(1).datafields);
set(handles.imtypepopup, 'String', {'dots', 'image'});
set(handles.pointtypepopup, 'String', {'marked', 'guessed'});
set(handles.clustertypepopup, 'String', [{'All'} {sc.name}]);
setupData;
handles = updateSettings(hObject, handles);

% Update handles structure
set(hObject, 'DeleteFcn', @cleanUpGlobalsCallback);
guidata(hObject, handles);

% UIWAIT makes sminspector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sminspector_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on selection change in xtypepopup.
function xtypepopup_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to xtypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xtypepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xtypepopup
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xtypepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xtypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ytypepopup.
function ytypepopup_Callback(hObject, eventdata, handles)
% hObject    handle to ytypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ytypepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ytypepopup
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ytypepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ytypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in imtypepopup.
function imtypepopup_Callback(hObject, eventdata, handles)
% hObject    handle to imtypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imtypepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imtypepopup
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function imtypepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imtypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in clustertypepopup.
function clustertypepopup_Callback(hObject, eventdata, handles)
% hObject    handle to clustertypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clustertypepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clustertypepopup
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function clustertypepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustertypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tracknumberpopup.
function tracknumberpopup_Callback(hObject, eventdata, handles)
% hObject    handle to tracknumberpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tracknumberpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tracknumberpopup
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tracknumberpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracknumberpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pointtypepopup.
function pointtypepopup_Callback(hObject, eventdata, handles)
% hObject    handle to pointtypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pointtypepopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pointtypepopup
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pointtypepopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pointtypepopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setupData ()
global SMI_SM;
global SMI_TRACKS;
global SMI_DATA;
global SMI_KNOWN;
global SMI_GUESSED;

disp ('loading model data and fitting paths');
[SMI_DATA, SMI_KNOWN] = SMI_SM.createHMMData(SMI_TRACKS);
SMI_GUESSED = cell(size(SMI_TRACKS));
for j = 1:length(SMI_TRACKS)
    SMI_GUESSED{j} = SMI_SM.doViterbi(SMI_TRACKS(j));
end
disp ('done');

function handles = updateSettings(hObject, handles)
handles.xtype = get(handles.xtypepopup, 'Value');
handles.ytype = get(handles.ytypepopup, 'Value');
handles.imtype = get(handles.imtypepopup, 'Value');
handles.clustertype = get(handles.clustertypepopup, 'Value') - 1;
handles.tracknumber = get(handles.tracknumberpopup, 'Value') - 1;
handles.pointtype = get(handles.pointtypepopup, 'Value');
handles.markersize = str2double(get(handles.markersizeedit, 'String'));
guidata(hObject, handles);
plotCluster(hObject, handles);

function plotCluster (hObject, handles)
cla (handles.imAxes);

global SMI_DATA;
global SMI_GUESSED;
global SMI_KNOWN;
global SMI_SM;

switch (handles.pointtype)
    case 1
        ptsrc = SMI_KNOWN;
    case 2
        ptsrc = SMI_GUESSED;
end

if (handles.tracknumber == 0) %all tracks
    data = cell2mat(SMI_DATA);
    ptlabel = cell2mat(ptsrc);
else
    data = SMI_DATA{handles.tracknumber};
    ptlabel = ptsrc{handles.tracknumber};
end

sc = SMI_SM.segmentationClusters;
c = 'bgrcmyk'; c = [c c c c c c c];
ms = handles.markersize;

if (handles.imtype == 1) %points
    if (handles.clustertype == 0) %allclusters
        valid = false(size(sc));
        for j = 1:length(sc)
            plot (handles.imAxes, data(handles.xtype, ptlabel == j), data (handles.ytype, ptlabel == j), [c(j) '.'], 'MarkerSize', ms);
            valid(j) = any(ptlabel == j);
            hold (handles.imAxes, 'on');
        end
        legend(handles.imAxes, {sc(valid).name}, 'Location', 'Best');
    else
        plot (handles.imAxes, data(handles.xtype, ptlabel > 0), data (handles.ytype, ptlabel > 0),'k.',...
            'Color', [0.7 0.7 0.7], 'MarkerSize', ms);
        hold (handles.imAxes, 'on');

        j = handles.clustertype;
        plot (handles.imAxes, data(handles.xtype, ptlabel == j), data (handles.ytype, ptlabel == j), [c(j) '.'], 'MarkerSize', ms);
        legend(handles.imAxes, {'All', sc(j).name}, 'Location', 'Best');        
    end
end
            
       
    

function cleanUpGlobalsCallback(src, eventdata)
clear global SMI_SM;
clear global SMI_TRACKS;
clear global SMI_DATA;
clear global SMI_KNOWN;
clear global SMI_GUESSED;



function markersizeedit_Callback(hObject, eventdata, handles)
% hObject    handle to markersizeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of markersizeedit as text
%        str2double(get(hObject,'String')) returns contents of markersizeedit as a double
updateSettings(hObject, handles);

% --- Executes during object creation, after setting all properties.
function markersizeedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markersizeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
