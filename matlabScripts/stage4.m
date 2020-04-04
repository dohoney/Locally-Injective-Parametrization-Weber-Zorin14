function varargout = stage4(varargin)
% STAGE4 MATLAB code for stage4.fig
%      STAGE4, by itself, creates a new STAGE4 or raises the existing
%      singleton*.
%
%      H = STAGE4 returns the handle to a new STAGE4 or the handle to
%      the existing singleton*.
%
%      STAGE4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STAGE4.M with the given input arguments.
%
%      STAGE4('Property','Value',...) creates a new STAGE4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stage4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stage4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stage4

% Last Modified by GUIDE v2.5 11-Oct-2015 15:22:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stage4_OpeningFcn, ...
                   'gui_OutputFcn',  @stage4_OutputFcn, ...
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


% --- Executes just before stage4 is made visible.
function stage4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stage4 (see VARARGIN)

% Choose default command line output for stage4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stage4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stage4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in showMesh.
function showMesh_Callback(hObject, eventdata, handles)
% hObject    handle to showMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
matrix = evalin('base','matrix');
tri_indices = evalin('base','tri_indices');
tri_faces=(reshape(tri_indices,3,length(tri_indices)/3))';
tri_faces = tri_faces+1;
x = real(matrix);
y = imag(matrix);
trimesh(tri_faces,x,y);
axis equal


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcbf)
