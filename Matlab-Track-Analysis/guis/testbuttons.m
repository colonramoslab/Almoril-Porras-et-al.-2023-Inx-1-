function varargout = testbuttons(varargin)
% TESTBUTTONS M-file for testbuttons.fig
%      TESTBUTTONS, by itself, creates a new TESTBUTTONS or raises the existing
%      singleton*.
%
%      H = TESTBUTTONS returns the handle to a new TESTBUTTONS or the handle to
%      the existing singleton*.
%
%      TESTBUTTONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTBUTTONS.M with the given input arguments.
%
%      TESTBUTTONS('Property','Value',...) creates a new TESTBUTTONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before testbuttons_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to testbuttons_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help testbuttons

% Last Modified by GUIDE v2.5 12-May-2010 15:36:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testbuttons_OpeningFcn, ...
                   'gui_OutputFcn',  @testbuttons_OutputFcn, ...
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


% --- Executes just before testbuttons is made visible.
function testbuttons_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testbuttons (see VARARGIN)

% Choose default command line output for testbuttons
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testbuttons wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testbuttons_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for j = 1:100
    set(handles.text1, 'String', num2str(j));
    if (get(hObject,'Value') == false)
        break;
    end
    pause(0.02);
end

% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1
for j = 1001:1100
    set(handles.text1, 'String',num2str(j));
    if (get(hObject,'Value') == false)
        break;
    end
    pause(0.02);
end
