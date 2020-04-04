function varargout = stage2(varargin)
% STAGE2 MATLAB code for stage2.fig
%      STAGE2, by itself, creates a new STAGE2 or raises the existing
%      singleton*.
%
%      H = STAGE2 returns the handle to a new STAGE2 or the handle to
%      the existing singleton*.
%
%      STAGE2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STAGE2.M with the given input arguments.
%
%      STAGE2('Property','Value',...) creates a new STAGE2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stage2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stage2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stage2

% Last Modified by GUIDE v2.5 10-Jan-2016 11:38:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stage2_OpeningFcn, ...
                   'gui_OutputFcn',  @stage2_OutputFcn, ...
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


% --- Executes just before stage2 is made visible.
function stage2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stage2 (see VARARGIN)

% Choose default command line output for stage2
handles.output = hObject;
handles.is3D = 1;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stage2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stage2_OutputFcn(hObject, eventdata, handles) 
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
m_faces = evalin('base','m_faces');
m_points = evalin('base','m_points');
%mn_points = (m_points-min(min(m_points)))/(max(max(m_points))-min(min(m_points))); 
 
trimesh(m_faces,m_points(:,1),m_points(:,2),m_points(:,3));


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(gcbf)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button2D.
function button2D_Callback(hObject, eventdata, handles)
% hObject    handle to button2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button2D
m_faces = evalin('base','m_faces');
m_points = evalin('base','m_points'); 
if handles.is3D==1
    trimesh(m_faces,m_points(:,1),m_points(:,2));
    handles.is3D = 0;
else
    trimesh(m_faces,m_points(:,1),m_points(:,2),m_points(:,3));
    handles.is3D = 1;
end
guidata(hObject, handles);

