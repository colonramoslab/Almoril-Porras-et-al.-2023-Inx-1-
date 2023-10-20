function varargout = TrackPlayer(varargin)
% TRACKPLAYER M-file for TrackPlayer.fig
%      TRACKPLAYER, by itself, creates a new TRACKPLAYER or raises the existing
%      singleton*.
%
%      H = TRACKPLAYER returns the handle to a new TRACKPLAYER or the handle to
%      the existing singleton*.
%
%      TRACKPLAYER('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in TRACKPLAYER.M with the given input arguments.
%
%      TRACKPLAYER('Property','Value',...) creates a new TRACKPLAYER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrackPlayer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrackPlayer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrackPlayer

% Last Modified by GUIDE v2.5 13-Jul-2010 11:38:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrackPlayer_OpeningFcn, ...
                   'gui_OutputFcn',  @TrackPlayer_OutputFcn, ...
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


% --- Executes just before TrackPlayer is made visible.
function TrackPlayer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrackPlayer (see VARARGIN)

% Choose default command line output for TrackPlayer
handles.output = hObject;
handles.ptbuffer = 50;
handles.imSize = 40;
handles.trackfield = 'sloc';

[handles,varargin] = assignToStruct(handles,varargin); 

global tp_alltracks;

if (isempty(varargin))
    error ('no args passed. TrackPlayer(src, options), where src is eset(s), expt(s), track(s)');
end
src = varargin{1};
if (isa (src, 'ExperimentSet'))
    e = [src.expt];
    tp_alltracks= [e.track];
else
    if (isa (src, 'Experiment'))
        tp_alltracks = [src.track];
    else
        if (isa (src, 'Track'))
            tp_alltracks = src;
        else
            error ('TrackPlayer(src, options), where src is eset(s), expt(s), track(s)');
        end
    end
end

handles.track = tp_alltracks(1);

set(handles.trackNumPopUp, 'String', cellfun(@num2str, num2cell(1:length(tp_alltracks)),'UniformOutput', false));


%load data type list into popup menu
topfields = {'speed', 'deltatheta', 'scovRatio'};

fieldList = [topfields, setdiff(unique([fieldnames(handles.track.dq); (handles.track.validDQName)']), topfields)];
handles.currentField = fieldList{1};
set(handles.dataType, 'String', fieldList);
handles = newTrack(handles, 1);
handles.lastButtonPressed = handles.play;
handles.busy = false;
% Update handles structure
set(handles.moviePlayButton, 'UserData', false);
set(handles.moviePlayButton, 'String', 'Play');
set(hObject, 'DeleteFcn', @cleanUpGlobalsCallback);
guidata(hObject, handles);

title (handles.imAxes, 'Let''s party!');

% UIWAIT makes TrackPlayer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrackPlayer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in dataType.
function dataType_Callback(hObject, eventdata, handles)
% hObject    handle to dataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dataType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataType

contents = cellstr(get(hObject,'String'));
handles.currentField = contents{get(hObject,'Value')};
handles = drawDataField(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dataType_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to dataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function startFrame_Callback(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrame as text
%        str2double(get(hObject,'String')) returns contents of startFrame as a double
handles.sf = str2double(get(hObject,'String'));
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function startFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endFrame_Callback(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrame as text
%        str2double(get(hObject,'String')) returns contents of endFrame as a double
handles.ef = str2double(get(hObject,'String'));
handles = rangeCheck(handles);
handles.cf = handles.ef;
handles = drawImage(handles);
handles = drawDataField(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function endFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in startDecrease.
function startDecrease_Callback(hObject, eventdata, handles)
% hObject    handle to startDecrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sf = handles.sf - 1;
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);


% --- Executes on button press in endDecrease.
function endDecrease_Callback(hObject, eventdata, handles)
% hObject    handle to endDecrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ef = handles.ef - 1;
handles = rangeCheck(handles);
handles.cf = handles.ef;
handles = drawImage(handles);
handles = drawDataField(handles);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);



% --- Executes on button press in startIncrease.
function startIncrease_Callback(hObject, eventdata, handles)
% hObject    handle to startIncrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sf = handles.sf + 1;
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);


% --- Executes on button press in endIncrease.
function endIncrease_Callback(hObject, eventdata, handles)
% hObject    handle to endIncrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ef = handles.ef + 1;
handles = rangeCheck(handles);
handles.cf = handles.ef;
handles = drawImage(handles);
handles = drawDataField(handles);
guidata(hObject, handles);

% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%movie should take at least 2 seconds to play
%{
delayTime = 2 / (handles.ef - handles.sf + 10);
delayTime = max(delayTime, 0.0625);
%}
delayTime = str2double(get(handles.delayTimeEdit, 'String'));
minDelay = 0.01;
ind0 = handles.interpedtopts(max(handles.sf - 5, 1));
ind1 = handles.interpedtopts(min(handles.ef + 5, length(handles.interpedtopts)));
inds = handles.ptstointerped(ind0:ind1);

for cf = inds
    ts = tic;
    handles.cf = cf;
    drawImage(handles);
    pause(max(delayTime - toc(ts), minDelay));
end
guidata(hObject, handles);

% --- Executes on button press in zoomOut.
%     resets pathAxes to maximum extent of track
function zoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to zoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loc = handles.track.getDerivedQuantity(handles.trackfield);
minl = min(loc,[],2);
maxl = max(loc,[],2);
xl = get(handles.pathAxes, 'XLim');
yl = get(handles.pathAxes, 'YLim');
xl = 1.5 * xl(1:2) - 0.5*xl(2:-1:1);
yl = 1.5 * yl(1:2) - 0.5*yl(2:-1:1);
xl = max(xl, minl(1));
xl = min(xl, maxl(1));
yl = max(yl, minl(2));
yl = min(yl, maxl(2));

axis(handles.pathAxes, [xl(1), xl(2), yl(1), yl(2)]);
axis(handles.pathAxes, 'equal');
axis(handles.pathAxes, [xl(1), xl(2), yl(1), yl(2)]);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);


% --- Executes on button press in zoomIn.
function zoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to zoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
r = getrect(handles.pathAxes);
axis(handles.pathAxes, [r(1) r(1)+r(3) r(2) r(2)+r(4)]);
axis(handles.pathAxes, 'equal');
axis(handles.pathAxes, [r(1) r(1)+r(3) r(2) r(2)+r(4)]);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);

% --- Executes on button press in pickPoints.
function pickPoints_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% hObject    handle to pickPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = getpts(handles.pathAxes);
[~,ind] = handles.track.nearestPoint([x,y]');
ind = handles.ptstointerped(ind);
%ind = handles.track.getDerivedQuantity('mapptstointerped', false, ind);
handles.sf = min(ind);
handles.ef = max(ind);
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
handles.lastButtonPressed = hObject;
handles = drawTrack(handles, true);
handles = tpfromcf(handles);
guidata(hObject, handles);



%plots the data field vs. time in dataAxes
function handles = drawDataField(handles, detailed)
existsAndDefault('detailed', true);
cla(handles.dataAxes);
s = max(handles.sf - handles.ptbuffer, 1);
e = min(handles.ef + handles.ptbuffer, length(handles.track.getDerivedQuantity('eti')));
handles.track.plotFields('eti', handles.currentField, 'inds', handles.sf:handles.ef, 'makeLegend', false, 'Axes', handles.dataAxes, 'LineWidth', 3, 'Color', 'y');
hold(handles.dataAxes, 'on');
handles.track.plotFields('eti', handles.currentField, 'inds', s:e, 'labeled', detailed, 'Axes', handles.dataAxes);

%draws the path in the pathAxes
function handles = drawTrack(handles, keepAxes)
existsAndDefault('keepAxes', false);
if (keepAxes)
    props = {'XLim', 'YLim'};
    vals = get(handles.pathAxes, props);
end

cla(handles.pathAxes);
handles.track.plotSegmentation('fieldName', handles.trackfield, 'multicolor', false, 'Axes', handles.pathAxes); hold(handles.pathAxes, 'on');
handles.track.plotPath(handles.trackfield, 'r-', 'inds', handles.sf:handles.ef, 'Axes', handles.pathAxes); hold off

axis(handles.pathAxes, 'equal');
hold (handles.pathAxes, 'off');
if (keepAxes)
    set(handles.pathAxes, props, vals);
end


%draws the image & overlays the track in imAxes
function handles = drawImage(handles)
cla(handles.imAxes);
ind = handles.interpedtopts(handles.cf);
pt = handles.pt;
pt(ind).drawTrackImage([],'Axes',handles.imAxes, 'fid', handles.track.expt.fid);
if (handles.cf < handles.sf || handles.cf > handles.ef)
    colormap gray;
else
    colormap (0.3*hot + 0.7*gray);
end
axis(handles.imAxes, 'equal');
yl = get(handles.imAxes, 'YLim');
xl = get(handles.imAxes, 'XLim');
xc = mean(xl);
yc = mean(yl);
if (diff(yl) > diff(xl))
    xw = handles.imSize;
    yw = xw * diff(yl) / diff(xl);
else
    yw = handles.imSize;
    xw = yw * diff(xl) / diff(yl);
end
xl = xc + xw * [-0.5 0.5];
yl = yc + yw * [-0.5 0.5];


hold (handles.imAxes, 'on');
s = max(handles.sf - handles.ptbuffer, 1);
e = min(handles.ef + handles.ptbuffer, length(handles.track.getDerivedQuantity('eti')));
handles.track.plotPath(handles.trackfield, 'b-', 'inds', s:e, 'highlightinds', handles.sf:handles.ef,...
    'highlightlinetype', 'r-', 'Axes', handles.imAxes);
hold (handles.imAxes, 'off');
set(handles.imAxes, 'XLim', xl, 'YLim', yl,'XTick',[],'YTick',[]);

if (handles.track.isrun(handles.cf))
    title (handles.imAxes, 'run');
else
    if (any(handles.track.sharpTurn.containsIndex(handles.cf)))
        title (handles.imAxes, handles.track.sharpTurn(handles.track.sharpTurn.containsIndex(handles.cf)).type);
    else
        if (any(handles.track.reorientation.containsIndex(handles.cf)))
            title (handles.imAxes, 'reorientation');
        else
            title (handles.imAxes, 'unknown');
        end
    end
end

function handles = rangeCheck(handles)
handles.sf = max(handles.sf, 1);
handles.ef = min(handles.ef, length(handles.track.getDerivedQuantity('eti')));
handles.ef = max(handles.ef, 1);
handles.sf = min(handles.sf, handles.ef);
set(handles.endFrame, 'String', num2str(handles.ef));
set(handles.startFrame, 'String', num2str(handles.sf));
set(handles.currentFrameEdit, 'String', num2str(handles.cf));


% --- Executes on button press in pickPointsDA.
function pickPointsDA_Callback(hObject, eventdata, handles)
% hObject    handle to pickPointsDA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = getpts(handles.dataAxes);
x = interp1(handles.track.getDerivedQuantity('eti'), 1:length(handles.track.getDerivedQuantity('eti')),x, 'nearest');
handles.sf = min(x);
handles.ef = max(x);
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
handles = drawTrack(handles, true);
handles.lastButtonPressed = hObject;
handles = tpfromcf(handles);
guidata(hObject, handles);



% --- Executes on selection change in trackNumPopUp.
function trackNumPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to trackNumPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trackNumPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trackNumPopUp
tracknum = get(hObject,'Value');
handles = newTrack(handles, tracknum);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function trackNumPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackNumPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = newTrack(handles, tracknum)
global tp_alltracks;
handles.track = tp_alltracks(tracknum);
handles.tracknum = tracknum;
handles.track.expt.openDataFile;
handles.ptstointerped = handles.track.getDerivedQuantity('mapptstointerped');
handles.interpedtopts = handles.track.getDerivedQuantity('mapinterpedtopts');
handles.pt = handles.track.pt;
handles.sf = 1;
handles.cf = 2;
handles.ef = 3;
fn = fieldnames(handles.track);
tpfields = {};
for j = 1:length(fn)
    if (~isempty(handles.track.(fn{j})) && isa (handles.track.(fn{j}), 'TrackPart'))
        tpfields = [tpfields, fn{j}]; %#ok<AGROW>
    end
end
set(handles.tptypePopUp, 'String', tpfields);
set(handles.tptypePopUp, 'Value', 1);
handles.tpn = 1;
handles = tpfromcf(handles);
handles = tpupdate(handles);
%set(handles.tpnumEdit, 'String', num2str(1));
%handles = rangeCheck(handles);
%handles = tprangeCheck(handles);
drawDataField(handles);
drawTrack(handles);
drawImage(handles);




% --- Executes on button press in fastforward.
function fastforward_Callback(hObject, eventdata, handles)
% hObject    handle to fastforward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fastforward
if (handles.busy)
    return;
end
delta = 1;

handles.busy = true;
handles = ffrwButtonAction(hObject, handles);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);

handles = continuousMovie(hObject, handles, delta);
handles.busy = false;

guidata(hObject, handles);


% --- Executes on button press in rewind.
function rewind_Callback(hObject, eventdata, handles)
% hObject    handle to rewind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rewind
if (handles.busy)
    return;
end
delta = -1;

handles.busy = true;
handles = ffrwButtonAction(hObject, handles);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);

handles = continuousMovie(hObject, handles, delta);
handles.busy = false;

guidata(hObject, handles);


function handles = continuousMovie (hObject, handles, delta)
%function continuousMovie (hObject, handles, delta)
%delayTime = .0625;
delayTime = str2double(get(handles.delayTimeEdit, 'String'));
minDelay = 0.01;
maxtime = 120;
tstart = tic;

while (get(hObject, 'Value') == true && toc(tstart) < maxtime)
   ts = tic;
   nextpt = min(handles.interpedtopts(handles.cf) + delta, handles.interpedtopts(end));
   nextpt = max(nextpt, handles.interpedtopts(1));
   handles.cf = handles.ptstointerped(nextpt);
   handles.(handles.edgefield) = handles.cf;
   handles.sf = min(handles.sf, handles.cf);
   handles.ef = max(handles.ef, handles.cf);
   handles = rangeCheck(handles);
   handles = drawImage(handles);
   handles = drawDataField(handles, false);
   pause(max(delayTime - toc(ts), minDelay));
end
handles = tpfromcf(handles);

function handles = ffrwButtonAction(hObject, handles)

switch (handles.lastButtonPressed)
    case handles.play
        if(hObject == handles.fastforward)
            handles.sf = handles.ef;
            handles.cf = handles.ef;
            handles.edgefield = 'ef';
        else
            handles.cf = handles.sf;
            handles.edgefield = 'sf';
        end
    case {handles.endIncrease, handles.endDecrease}
        handles.cf = handles.ef;
        handles.edgefield = 'ef';
    case {handles.startIncrease, handles.startDecrease}
        handles.cf = handles.sf;
        handles.edgefield = 'sf';
    case handles.goforward
        handles.edgefield = 'ef';
    case {handles.fastforward, handles.rewind}
        %don't touch
    otherwise
        if (hObject == handles.fastforward)
            handles.cf = handles.ef;
            handles.edgefield = 'ef';
        else
            handles.cf = handles.sf;
            handles.edgefield = 'sf';
        end
end


% --- Executes on button press in pickPointsIm.
function pickPointsIm_Callback(hObject, eventdata, handles)
% hObject    handle to pickPointsIm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,y] = getpts(handles.imAxes);
[~,ind] = handles.track.nearestPoint([x y]');
handles.sf = min(handles.ptstointerped(ind));
handles.ef = max(handles.ptstointerped(ind));
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
handles.lastButtonPressed = hObject;
handles = tpfromcf(handles);
guidata(hObject, handles);


function cleanUpGlobalsCallback(src, eventdata)
clear global tp_alltracks;


function delayTimeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to delayTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delayTimeEdit as text
%        str2double(get(hObject,'String')) returns contents of delayTimeEdit as a double


% --- Executes during object creation, after setting all properties.
function delayTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delayTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in moviePlayButton.
function moviePlayButton_Callback(hObject, eventdata, handles)
% hObject    handle to moviePlayButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~isempty(get(hObject, 'UserData')) && get(hObject, 'UserData'))
    set(hObject, 'UserData', false);
    set(hObject, 'String', 'Play');
    return
end

set(hObject, 'UserData', true);
set(hObject, 'String', 'Pause');
df = handles.ef - handles.cf;
db = handles.cf - handles.sf;
cla(handles.dataAxes);
handles.track.plotFields('eti', handles.currentField, 'labeled', true, 'Axes', handles.dataAxes); hold(handles.dataAxes, 'on');
%h = [];
%eti = handles.track.getDerivedQuantity('eti');
%ydat = handles.track.getDerivedQuantity(handles.currentField);
update = 0;
oldval = get(handles.dataType, 'Value');
ts = tic;
while (get(hObject, 'UserData'))

    if (oldval ~= get(handles.dataType, 'Value'))
        update = 2;
        contents = cellstr(get(handles.dataType,'String'));
        handles.currentField = contents{get(handles.dataType,'Value')};
        oldval = get(handles.dataType, 'Value');
    else
        update = update - 1;
    end
    delayTime = str2double(get(handles.delayTimeEdit, 'String'));
    handles.cf = handles.cf + 1;
    handles.sf = handles.cf - db;
    handles.ef = handles.cf + df;
    
    handles = rangeCheck(handles);
    handles = drawImage(handles);
   %{
    if (~isempty(h))
        try
            delete(h);
        catch %#ok<CTCH>
        end
    end
    %}
   % h = handles.track.plotFields('eti', handles.currentField, 'inds', handles.sf:handles.ef, 'makeLegend', false, 'Axes', handles.dataAxes, 'LineWidth', 3, 'Color', 'y'); hold on
  %  h = plot (handles.dataAxes,eti(handles.cf), ydat(:,handles.cf), 'co', 'MarkerSize', 10, 'LineWidth', 2);
    xlim(handles.dataAxes, handles.track.getDerivedQuantity('eti', false, [handles.sf handles.ef]));
    pause(max(delayTime - toc(ts), 0.001));
    set(handles.totalTime, 'String', num2str(toc(ts), 3));

    if (update > 0)
        cla(handles.dataAxes);
        handles.track.plotFields('eti', handles.currentField, 'labeled', true, 'Axes', handles.dataAxes); hold(handles.dataAxes, 'on');
     %   h = [];
      %  ydat = handles.track.getDerivedQuantity(handles.currentField);
       
    end
    ts = tic;

end
hold (handles.dataAxes, 'off');
handles = tpfromcf(handles);
guidata(hObject, handles);



function currentFrameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to currentFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentFrameEdit as text
%        str2double(get(hObject,'String')) returns contents of currentFrameEdit as a double
handles.cf = str2double(get(hObject,'String'));
handles = rangeCheck(handles);
handles = drawImage(handles);
handles = drawDataField(handles);
handles = tpfromcf(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function currentFrameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reportEdit_Callback(hObject, eventdata, handles)
% hObject    handle to reportEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reportEdit as text
%        str2double(get(hObject,'String')) returns contents of reportEdit as a double


% --- Executes during object creation, after setting all properties.
function reportEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reportEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in tptypePopUp.
function tptypePopUp_Callback(hObject, eventdata, handles)
% hObject    handle to tptypePopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tptypePopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tptypePopUp
handles = tpfromcf(handles);
handles = tpupdate(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tptypePopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tptypePopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tpnumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to tpnumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tpnumEdit as text
%        str2double(get(hObject,'String')) returns contents of tpnumEdit as a double
handles.tpn = str2double(get(hObject,'String'));
handles = tpupdate(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tpnumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tpnumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tpnextButton.
function tpnextButton_Callback(hObject, eventdata, handles)
% hObject    handle to tpnextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tpn = handles.tpn + 1;
handles = tpupdate(handles);
guidata(hObject, handles);


% --- Executes on button press in tpprevButton.
function tpprevButton_Callback(hObject, eventdata, handles)
% hObject    handle to tpprevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.tpn = handles.tpn - 1;
handles = tpupdate(handles);
guidata(hObject, handles);


% --- Executes on button press in segreportButton.
function segreportButton_Callback(hObject, eventdata, handles)
% hObject    handle to segreportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.reportEdit, 'String', handles.track.getReport(handles.sf, handles.ef));


function handles = tprangeCheck(handles)
fn = get(handles.tptypePopUp, 'String');
fn = fn{get(handles.tptypePopUp, 'Value')};
tpn = handles.tpn;
tpl = length(handles.track.(fn));
set(handles.tpmaxString, 'String', num2str(tpl));
tpn = max(tpn, 1);
tpn = min(tpn, tpl);
set(handles.tpnumEdit, 'String', num2str(tpn));
handles.tpn = tpn;

function handles = tpupdate(handles)
handles = tprangeCheck(handles);
fn = get(handles.tptypePopUp, 'String');
fn = fn{get(handles.tptypePopUp, 'Value')};
tp = handles.track.(fn)(handles.tpn);
handles.sf = tp.startInd;
handles.ef = tp.endInd;
if (isfield(tp, 'centralInd'))
    handles.cf = round(tp.centralInd);
else
    handles.cf = round(mean([handles.sf handles.ef]));
end
handles = rangeCheck(handles);
handles = drawImage(handles);
handles = drawDataField(handles);
set(handles.reportEdit, 'String', tp.getReport);

function handles = tpfromcf (handles)
fn = get(handles.tptypePopUp, 'String');
fn = fn{get(handles.tptypePopUp, 'Value')};
tpn = find(handles.track.(fn).containsIndex(handles.cf), 1, 'first');
if (~isempty(tpn))
    handles.tpn = tpn;
else
    tpn = find([handles.track.(fn).startInd] >= handles.sf, 1, 'first');
    if (~isempty(tpn))
        handles.tpn = tpn;
    else
        handles.tpn = length(handles.(fn));
    end
end
handles = tprangeCheck(handles);
