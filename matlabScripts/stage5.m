function varargout = stage5(varargin)
% STAGE5 MATLAB code for stage5.fig
%      STAGE5, by itself, creates a new STAGE5 or raises the existing
%      singleton*.
%
%      H = STAGE5 returns the handle to a new STAGE5 or the handle to
%      the existing singleton*.
%
%      STAGE5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STAGE5.M with the given input arguments.
%
%      STAGE5('Property','Value',...) creates a new STAGE5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stage5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stage5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stage5

% Last Modified by GUIDE v2.5 14-Jan-2016 14:27:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stage5_OpeningFcn, ...
                   'gui_OutputFcn',  @stage5_OutputFcn, ...
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


% --- Executes just before stage5 is made visible.
function stage5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stage5 (see VARARGIN)

% Choose default command line output for stage5
handles.output = hObject;

axes(handles.axes5);
imshow('arrow.png');
axes(handles.axes6);
imshow('arrow.png');

handles.render = 0;

% Update handles structure
guidata(hObject, handles);
%axes(handles.axes3);

% UIWAIT makes stage5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stage5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.render == 1
    delete(gcbf)
end

uSource = evalin('base','uSource');
weightsMatSource = evalin('base','weightsMatSource');
s_t=tic;
outSource = weightsMatSource\uSource;
sTime = toc(s_t);

uTarget = evalin('base','uTarget');
weightsMatTarget = evalin('base','weightsMatTarget');
t_t = tic;
outTarget = weightsMatTarget\uTarget;
tTime = toc(t_t);

outSource = full(outSource);
outTarget = full(outTarget);
%b1=zeros(length(outSource),2);
%b2=zeros(length(outTarget),2);
%for i=1:length(outSource)
%    b1(i,1) = outSource(i,1);
%    b1(i,2) = outSource(i,2);
%end
%for i=1:length(outTarget)
%    b2(i,1) = outTarget(i,1);
%    b2(i,2) = outTarget(i,2);
%end
assignin('base','outSource',outSource);
assignin('base','outTarget',outTarget);
assignin('base','sTime',sTime);
assignin('base','tTime',tTime);

delete(gcbf)


% --- Executes on button press in renderButton.
function renderButton_Callback(hObject, eventdata, handles)
% hObject    handle to renderButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

m_points = evalin('base','m_points');
m_faces = evalin('base','m_faces');
axes(handles.axes1);
trimesh(m_faces,m_points(:,1),m_points(:,2),m_points(:,3));
axes(handles.axes2);
matrix = evalin('base','matrix');
tri_indices = evalin('base','tri_indices');
tri_faces=(reshape(tri_indices,3,length(tri_indices)/3))';
tri_faces = tri_faces+1;
x = real(matrix);
y = imag(matrix);
trimesh(tri_faces,x,y);
axis equal
axes(handles.axes3);
uSource = evalin('base','uSource');
weightsMatSource = evalin('base','weightsMatSource');
s_t=tic;
outSource = weightsMatSource\uSource;
sTime = toc(s_t);
trimesh(m_faces,outSource(:,1),outSource(:,2))
axis equal
axes(handles.axes4);
uTarget = evalin('base','uTarget');
weightsMatTarget = evalin('base','weightsMatTarget');
t_t = tic;
outTarget = weightsMatTarget\uTarget;
tTime = toc(t_t);
trimesh(tri_faces,outTarget(:,1),outTarget(:,2))
axis equal
b1=zeros(length(outSource),2);
b2=zeros(length(outTarget),2);
for i=1:length(outSource)
    b1(i,1) = outSource(i,1);
    b1(i,2) = outSource(i,2);
end
for i=1:length(outTarget)
    b2(i,1) = outTarget(i,1);
    b2(i,2) = outTarget(i,2);
end
assignin('base','outSource',b1);
assignin('base','outTarget',b2);
assignin('base','sTime',sTime);
assignin('base','tTime',tTime);

handles.render = 1;
guidata(hObject, handles);
