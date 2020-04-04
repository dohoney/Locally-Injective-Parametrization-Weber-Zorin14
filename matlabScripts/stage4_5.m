function varargout = stage4_5(varargin)
% STAGE4_5 MATLAB code for stage4_5.fig
%      STAGE4_5, by itself, creates a new STAGE4_5 or raises the existing
%      singleton*.
%
%      H = STAGE4_5 returns the handle to a new STAGE4_5 or the handle to
%      the existing singleton*.
%
%      STAGE4_5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STAGE4_5.M with the given input arguments.
%
%      STAGE4_5('Property','Value',...) creates a new STAGE4_5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stage4_5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stage4_5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stage4_5

% Last Modified by GUIDE v2.5 24-Jan-2016 20:54:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stage4_5_OpeningFcn, ...
                   'gui_OutputFcn',  @stage4_5_OutputFcn, ...
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


% --- Executes just before stage4_5 is made visible.
function stage4_5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stage4_5 (see VARARGIN)

% Choose default command line output for stage4_5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stage4_5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stage4_5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CC.
function CC_Callback(hObject, eventdata, handles)
% hObject    handle to CC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
weightsSelect = [1 1];
assignin('base','weightsSelect',weightsSelect);
delete(gcbf)

% --- Executes on button press in CM.
function CM_Callback(hObject, eventdata, handles)
% hObject    handle to CM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
weightsSelect = [1 0];
assignin('base','weightsSelect',weightsSelect);
delete(gcbf)

% --- Executes on button press in MC.
function MC_Callback(hObject, eventdata, handles)
% hObject    handle to MC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
weightsSelect = [0 1];
assignin('base','weightsSelect',weightsSelect);
delete(gcbf)

% --- Executes on button press in MM.
function MM_Callback(hObject, eventdata, handles)
% hObject    handle to MM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
weightsSelect = [0 0];
assignin('base','weightsSelect',weightsSelect);
delete(gcbf)