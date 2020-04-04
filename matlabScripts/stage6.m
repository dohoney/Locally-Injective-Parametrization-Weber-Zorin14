function varargout = stage6(varargin)
% STAGE6 MATLAB code for stage6.fig
%      STAGE6, by itself, creates a new STAGE6 or raises the existing
%      singleton*.
%
%      H = STAGE6 returns the handle to a new STAGE6 or the handle to
%      the existing singleton*.
%
%      STAGE6('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STAGE6.M with the given input arguments.
%
%      STAGE6('Property','Value',...) creates a new STAGE6 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stage6_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stage6_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stage6

% Last Modified by GUIDE v2.5 10-Jan-2016 12:46:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stage6_OpeningFcn, ...
                   'gui_OutputFcn',  @stage6_OutputFcn, ...
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


% --- Executes just before stage6 is made visible.
function stage6_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stage6 (see VARARGIN)

% Choose default command line output for stage6
handles.output = hObject;
handles.is3d = 1;
finalFvec = evalin('base','finalFvec');
finalOut = evalin('base','finalOut');
finalPvec = evalin('base','finalPvec');
axes(handles.sourceAxes);
trimesh( finalFvec , finalPvec(:,1) , finalPvec(:,2) , finalPvec(:,3) )
axes(handles.targetAxes);
trimesh( finalFvec , finalOut(:,1) , finalOut(:,2) )
axis equal


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stage6 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stage6_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in viewButton.
function viewButton_Callback(hObject, eventdata, handles)
% hObject    handle to viewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
m_faces = evalin('base','m_faces');
finalOut = evalin('base','finalOut');
trimesh( m_faces , finalOut(:,1) , finalOut(:,2) )
 axis 'equal'


% --- Executes on button press in exportButton.
function exportButton_Callback(hObject, eventdata, handles)
% hObject    handle to exportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
finalPvec = evalin('base','finalPvec');
finalFvec = evalin('base','finalFvec');
finalOut = evalin('base','finalOut');
obj.vertices = finalPvec;
obj.vertices_texture = finalOut;
objects.data.vertices = finalFvec;
objects.data.texture = finalFvec;
obj.objects = objects;
obj.objects.type='f';

FilterSpec = '*.obj';
[FileName,PathName] = uiputfile(FilterSpec);
if FileName == 0
    return;
end
write_wobj(obj,FileName);
cur = pwd;
cur =[cur '\'];
if strcmp(cur,PathName) ~= 1
    movefile(FileName,PathName,'f');
end

% --- Executes on button press in ex2222.
function ex2222_Callback(hObject, eventdata, handles)
% hObject    handle to ex2222 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m_points = evalin('base','m_points');
m_faces = evalin('base','m_faces');
finalOut = evalin('base','finalOut');
obj.vertices = m_points;
obj.vertices_texture = finalOut;
objects.data.vertices = m_faces;
objects.data.texture = m_faces;
obj.objects = objects;
obj.objects.type='f';

FilterSpec = '*.obj';
[FileName,PathName] = uiputfile(FilterSpec);
write_wobj(obj,FileName);
cur = pwd;
cur =[cur '\'];
if strcmp(cur,PathName) ~= 1
    movefile(FileName,PathName,'f');
end



% --- Executes on button press in button2D.
function button2D_Callback(hObject, eventdata, handles)
% hObject    handle to button2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button2D
handles.output = hObject;
finalFvec = evalin('base','finalFvec');
finalPvec = evalin('base','finalPvec');
axes(handles.sourceAxes);
if handles.is3d == 1
    trimesh( finalFvec , finalPvec(:,1) , finalPvec(:,2)  )
    handles.is3d = 0;
else
    trimesh( finalFvec , finalPvec(:,1) , finalPvec(:,2) , finalPvec(:,3) )
    handles.is3d = 1;
end

guidata(hObject, handles);
