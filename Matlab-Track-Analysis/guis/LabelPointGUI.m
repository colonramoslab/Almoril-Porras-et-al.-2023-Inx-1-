function varargout = LabelPointGUI(varargin)
% LABELPOINTGUI M-file for LabelPointGUI.fig
%      LABELPOINTGUI, by itself, creates a new LABELPOINTGUI or raises the existing
%      singleton*.
%
%      H = LABELPOINTGUI returns the handle to a new LABELPOINTGUI or the handle to
%      the existing singleton*.
%
%      LABELPOINTGUI('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in LABELPOINTGUI.M with the given input arguments.
%
%      LABELPOINTGUI('Property','Value',...) creates a new LABELPOINTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LabelPointGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LabelPointGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LabelPointGUI

% Last Modified by GUIDE v2.5 07-Jul-2010 11:16:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LabelPointGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LabelPointGUI_OutputFcn, ...
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


% --- Executes just before LabelPointGUI is made visible.
function LabelPointGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LabelPointGUI (see VARARGIN)

% Choose default command line output for LabelPointGUI
handles.output = hObject;
handles.ptbuffer = 50;
handles.trackfield = 'sloc';
sm = [];

%find the path where we stick saved segmentation clusters
mf = mfilename('fullpath');
thisbasedir = fileparts(fileparts(mf)); %go back behind "guis" folder
handles.basedir = fullfile(thisbasedir, 'Segmentation Models', '');
handles.username = [];

varargin = assignApplicable(varargin);
[handles,varargin] = assignToStruct(handles,varargin); %#ok<NASGU>

handles.markedtracks = [];
if (isdir(handles.basedir))
    if (handles.basedir(end) ~= filesep)
        handles.basedir = [handles.basedir filesep];
    end
    handles.lastsaved = [handles.basedir '*.mat'];
else
    handles.lastsaved = handles.basedir;
end
%system('git status')
if (isempty(sm) || ~isa (sm, 'SegmentationModel'))
    error ('you must pass a SegmentationModel object with ''sm'' flag');
end

global SEGMENTATION_MODEL;
global alltracks;

SEGMENTATION_MODEL = sm;

handles.sc = sm.segmentationClusters;
t = [sm.eset.expt.track];
[~,I] = sort([t.npts],'descend');
alltracks = t(I);
handles.track = alltracks(1);
set(handles.trackNumPopUp, 'String', cellfun(@num2str, num2cell(1:length(alltracks)),'UniformOutput', false));


%load data type list into popup menu
fieldList = [unique([handles.sc.datafields]) Track.validDQName];
handles.currentField = fieldList{1};
set(handles.dataType, 'String', fieldList);
set(handles.clusterName, 'String', {handles.sc.name});
set(handles.pointTypePopup, 'String', {'guessed points', 'marked points'}, 'Value', 2);
handles.pointtype = 'marked';
handles = newTrack(handles, 1);
handles.lastButtonPressed = handles.play;
handles.busy = false;
% Update handles structure
set(hObject, 'DeleteFcn', @cleanUpGlobalsCallback);
guidata(hObject, handles);

title (handles.imAxes, 'Let''s party!');

% UIWAIT makes LabelPointGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LabelPointGUI_OutputFcn(hObject, eventdata, handles) 
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
guidata(hObject, handles);


% --- Executes on selection change in clusterName.
function clusterName_Callback(hObject, eventdata, handles)
% hObject    handle to clusterName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clusterName contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusterName


% --- Executes during object creation, after setting all properties.
function clusterName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in markAs.
function markAs_Callback(hObject, eventdata, handles)
% hObject    handle to markAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ci = get(handles.clusterName,'Value');
answer = yesnodialog('title','Confirm Labeling Points', 'string', ['mark as ', handles.sc(ci).name]);
if (strcmpi(answer, 'yes'))
    handles.sc(ci).addKnownPoints(handles.track, handles.sf:handles.ef);
    set(handles.pointTypePopup, 'Value', 2);
    handles.pointtype = 'marked';
    props = {'XLim', 'YLim'};
    vals = get(handles.pathAxes, props);
    handles = drawTrack(handles);
    set(handles.pathAxes, props, vals);
    handles.markedtracks = unique([handles.markedtracks handles.tracknum]);
    
end
handles.lastButtonPressed = hObject;
guidata(hObject, handles);


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ci = get(handles.clusterName,'Value');
[~,mostlikely] = max(handles.sc.posteriorProbabilities(handles.track));
sf = find(1:length(mostlikely) > handles.ef & mostlikely == ci, 1, 'first');
if (~isempty(sf))
    handles.sf = sf;
end
ef = find(1:length(mostlikely) > handles.sf & mostlikely ~= ci, 1, 'first');
if (~isempty(ef))
    handles.ef = ef;
end
handles = rangeCheck(handles);
handles.cf = handles.sf;
handles = drawImage(handles);
handles = drawDataField(handles);
handles.lastButtonPressed = hObject;

guidata(hObject, handles);



%plots the data field vs. time in dataAxes
function handles = drawDataField(handles, detailed)
existsAndDefault('detailed', true);
cla(handles.dataAxes);
s = max(handles.sf - handles.ptbuffer, 1);
e = min(handles.ef + handles.ptbuffer, length(handles.track.getDerivedQuantity('eti')));
x = handles.track.getDerivedQuantity('eti', false, s:e);
y = handles.track.getDerivedQuantity(handles.currentField, false, s:e);
xh = handles.track.getDerivedQuantity('eti', false, handles.sf:handles.ef);
yh = handles.track.getDerivedQuantity(handles.currentField, false, handles.sf:handles.ef);
plot (handles.dataAxes,x,y,'k-',xh,yh,'r.-'); 
if (detailed)
    hold(handles.dataAxes, 'on');
    c = 'rgcmyk'; c = [c c c c c];
    msym = '******vvvvvv++++++xxxxxxhhhhhh';
    valid = false(size(handles.sc));
    for j = 1:length(handles.sc)
        if (strcmpi(handles.sc(j).name, 'run') || strcmpi(handles.sc(j).name, 'runs'))
            marker = 'b.';
        else
            marker = [c(j) msym(j)];
        end
        if (~isempty(handles.sc(j).knownPoints))
            [~,inds] = intersect(s:e,[handles.sc(j).knownPoints([handles.sc(j).knownPoints.track] == handles.track).inds]);
            if (~isempty(inds))
                plot (handles.dataAxes, x(inds), y(inds), marker);
                valid(j) = true;
            end
        end
    end
    legend(handles.dataAxes, ['path', 'current', {handles.sc(valid).name}], 'Location', 'Best');
    hold(handles.dataAxes, 'off');
end

%draws the path in the pathAxes
function handles = drawTrack(handles, keepAxes)
existsAndDefault('keepAxes', false);
if (keepAxes)
    props = {'XLim', 'YLim'};
    vals = get(handles.pathAxes, props);
end

cla(handles.pathAxes);
%handles.track.plotPath(handles.trackfield, 'b-', 'highlightinds', handles.sf:handles.ef, 'Axes', handles.pathAxes);
%axis(handles.pathAxes, 'equal');
ind = find(strcmpi({handles.sc.name}, 'run'));
if (~isempty(ind))
    indorder = ([ind setdiff(1:length(handles.sc), ind)]);
else
    ind = find(strcmpi({handles.sc.name}, 'runs'));
    if (~isempty(ind))
        indorder = ([ind setdiff(1:length(handles.sc), ind)]);
    else
        indorder = 1:length(handles.sc);
    end
end

c = 'rgcmyk'; c = [c c c c c];
msym = '******vvvvvv++++++xxxxxxhhhhhh';
switch (lower(handles.pointtype))
    case 'guessed'
        [~,mostlikely] = max(handles.sc.posteriorProbabilities(handles.track));
        handles.track.plotPath('sloc', 'b-', 'inds', mostlikely == indorder(1),'Axes', handles.pathAxes); 
        hold(handles.pathAxes, 'on');
        valid = false(size(indorder));
        for j = indorder(2:end);
            if (any(mostlikely == j))
                handles.track.plotPath('sloc', [c(j) '.'], 'inds', mostlikely == j,'Axes', handles.pathAxes);
                hold(handles.pathAxes, 'on');
                valid(j) = true;
            end
        end
        legend(handles.pathAxes, {handles.sc(indorder(valid)).name}, 'Location', 'Best');
    case 'marked'
        handles.track.plotPath('sloc', 'k-', 'Color', [0.5 0.5 0.5], 'Axes', handles.pathAxes);
        hold (handles.pathAxes, 'on');
        valid = false(size(indorder));
        for j = indorder
            if j == indorder(1)
                marker = 'b.';
            else
                marker = [c(j) msym(j)];
            end
            if (~isempty(handles.sc(j).knownPoints))
                inds = [handles.sc(j).knownPoints([handles.sc(j).knownPoints.track] == handles.track).inds];
                if (~isempty(inds))
                    handles.track.plotPath('sloc', marker, 'inds', inds,'Axes', handles.pathAxes);
                    hold(handles.pathAxes, 'on');
                    valid(j) = true;
                end            
            end
        end
        
        legh = legend(handles.pathAxes, ['path', {handles.sc(indorder(valid)).name}], 'Location', 'Best');
        
end
handles.track.plotPath('sloc', 'r-', 'inds', handles.sf:handles.ef, 'Axes', handles.pathAxes);

axis(handles.pathAxes, 'equal');
hold (handles.pathAxes, 'off');
if (keepAxes)
    set(handles.pathAxes, props, vals);
    set(legh, 'Location', 'BestOutside');
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
hold (handles.imAxes, 'on');
s = max(handles.sf - handles.ptbuffer, 1);
e = min(handles.ef + handles.ptbuffer, length(handles.track.getDerivedQuantity('eti')));
handles.track.plotPath(handles.trackfield, 'b-', 'inds', s:e, 'highlightinds', handles.sf:handles.ef,...
    'highlightlinetype', 'r-', 'Axes', handles.imAxes);
hold (handles.imAxes, 'off');
set(handles.imAxes, 'XLim', xl, 'YLim', yl,'XTick',[],'YTick',[]);
if (handles.sctype(handles.cf) > 0)
    title (handles.imAxes, handles.sc(handles.sctype(handles.cf)).name);
else
    title (handles.imAxes, '');
end

function handles = rangeCheck(handles)
handles.sf = max(handles.sf, 1);
handles.ef = min(handles.ef, length(handles.track.getDerivedQuantity('eti')));
handles.ef = max(handles.ef, 1);
handles.sf = min(handles.sf, handles.ef);
set(handles.endFrame, 'String', num2str(handles.ef));
set(handles.startFrame, 'String', num2str(handles.sf));


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
guidata(hObject, handles);


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if handles.basedir(end) ~= filesep
%    handles.basedir = [handles.basedir filesep];
%end
global SEGMENTATION_MODEL;
[fname, d] = uiputfile(handles.lastsaved, 'Enter a file name to save model');
if (fname ~= 0)
    SEGMENTATION_MODEL.toMatFile(fullfile(d,fname));
end
handles.lastsaved = fullfile(d,fname);
handles.lastButtonPressed = hObject;
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
global alltracks;
handles.track = alltracks(tracknum);
handles.tracknum = tracknum;
handles.track.expt.openDataFile;
handles.ptstointerped = handles.track.getDerivedQuantity('mapptstointerped');
handles.interpedtopts = handles.track.getDerivedQuantity('mapinterpedtopts');
handles.pt = handles.track.pt;
handles.sf = 1;
handles.ef = 2;
handles.cf = handles.sf;
handles = rangeCheck(handles);
handles = setSegmentationSequence(handles);
drawDataField(handles);
drawTrack(handles);
drawImage(handles);


% --- Executes on selection change in pointTypePopup.
function pointTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to pointTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pointTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pointTypePopup
switch get(hObject,'Value') 
    case 1
        handles.pointtype = 'guessed';
    case 2
        handles.pointtype = 'marked';
    otherwise
        handles.pointtype = 'marked';
end
props = {'XLim', 'YLim'};
vals = get(handles.pathAxes, props);
handles = drawTrack(handles);
set(handles.pathAxes, props, vals);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pointTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pointTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
maxtime = 30;
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
guidata(hObject, handles);


% --- Executes on button press in commitToGit.
function commitToGit_Callback(hObject, eventdata, handles)
% hObject    handle to commitToGit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.lastButtonPressed ~= handles.saveButton)
    answer = yesnodialog('title','Possible Unsaved Changes', 'string', 'unsaved changes: continue?');
    if (strcmpi(answer,'no'))
        return;
    end
end
username = getName('title', 'Enter Your Name For Commit', 'string', 'OK to commit; cancel to abort', 'default', handles.username);
if (isempty(username))
    return;
end
handles.username = username;
commitmsg = sprintf([username ': marked ' num2str(length(handles.markedtracks)) ' tracks']);
[pathstr, name, ext] = fileparts(handles.lastsaved);
pathstr = relativepath(pathstr);
fname = fullfile(pathstr, [name ext]);
%fname = handles.lastsaved;
fname(fname == '\') = '/';
%fname
%if (~isempty(ind))
%    fname(ind) = fname(ind-1);
%    fname(ind-1) = '/';
%end
command1 = ['git add "' fname '"'];
command2 = ['git commit -m "' commitmsg '"'];
[status, msg] = system(command1);
disp(msg);
if (status ~= 0)
    yesnodialog ('title', 'Error Adding File', 'string', 'see message in matlab output window');
else
    [status, msg] = system(command2);
    disp(msg);
    if (status ~= 0)
        yesnodialog ('title', 'Error Committing', 'string', 'see message in matlab output window');
    end
end
handles.markedtracks = [];


% --- Executes on button press in unmark.
function unmark_Callback(hObject, eventdata, handles)
% hObject    handle to unmark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for j = 1:length(handles.sc)
    valid = true(size(handles.sc(j).knownPoints));
    if (isempty(handles.sc(j).knownPoints))
        continue;
    end
    for k = find([handles.sc(j).knownPoints.track] == handles.track)
        handles.sc(j).knownPoints(k).inds = setdiff(handles.sc(j).knownPoints(k).inds, handles.sf:handles.ef);
        if (isempty(handles.sc(j).knownPoints(k).inds))
            valid(k) = false;
        end
    end
    handles.sc(j).knownPoints = handles.sc(j).knownPoints(valid);
end
handles = drawDataField(handles);
handles = drawTrack(handles, true);
handles.lastButtonPressed = hObject;
guidata(hObject, handles);


% --- Executes on button press in goforward.
function goforward_Callback(hObject, eventdata, handles)
% hObject    handle to goforward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SEGMENTATION_MODEL;
sm = SEGMENTATION_MODEL;
handles.sf = handles.ef + 1;
handles.ef = handles.sf + 1;
handles.cf = handles.ef;
handles = rangeCheck(handles);
handles = drawImage(handles);
handles = drawDataField(handles);
handles = drawTrack(handles, true);
handles.lastButtonPressed = hObject;
ci = get(handles.clusterName,'Value');
set(handles.clusterName, 'Value', sm.nextCluster(ci));
guidata(hObject, handles);

function cleanUpGlobalsCallback(src, eventdata)
clear global SEGMENTATION_MODEL;
clear global alltracks;

function handles = setSegmentationSequence(handles)
sctype = repmat(-1,size(handles.track.getDerivedQuantity('eti')));
for j = 1:length(handles.sc)
    kp = [handles.sc(j).knownPoints];
    if (~isempty(kp))
        kp = kp([kp.track] == handles.track);
    end
    if (~isempty(kp))
        sctype([kp.inds]) = j;
    end
end
handles.sctype = sctype;



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
