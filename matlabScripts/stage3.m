function varargout = stage3(varargin)
% STAGE3 MATLAB code for stage3.fig
%      STAGE3, by itself, creates a new STAGE3 or raises the existing
%      singleton*.
%
%      H = STAGE3 returns the handle to a new STAGE3 or the handle to
%      the existing singleton*.
%
%      STAGE3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STAGE3.M with the given input arguments.
%
%      STAGE3('Property','Value',...) creates a new STAGE3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stage3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stage3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stage3

% Last Modified by GUIDE v2.5 01-Feb-2016 13:00:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stage3_OpeningFcn, ...
                   'gui_OutputFcn',  @stage3_OutputFcn, ...
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


% --- Executes just before stage3 is made visible.
function stage3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stage3 (see VARARGIN)

% Choose default command line output for stage3
handles.output = hObject;
handles.new_border = [];
handles.toFlip = 0;
handles.edit = 0;
%handles.m_border = evalin('base','m_border');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stage3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stage3_OutputFcn(hObject, eventdata, handles) 
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

n_b = handles.bp;
n_bSize = length(n_b);
rotIndices = handles.rotIndices;
if handles.toFlip == 1
    n_b = flipud(n_b);
    rotIndices = flipud(rotIndices);
end
assignin('base','n_b',n_b);
assignin('base','n_bSize',n_bSize);
assignin('base','rotIndices',rotIndices);
delete(gcbf)


% --- Executes on button press in gen.
function gen_Callback(hObject, eventdata, handles)
% hObject    handle to gen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%O = evalin('base','O');
%bp = evalin('base','bp');
%m_border = evalin('base','m_border');
%m_border = flipud (m_border);

%m_border = (m_border-min(min(m_border)))/(max(max(m_border))-min(min(m_border)));
%plot ([m_points(O(:,1),1) m_points(O(:,2),1)]',[m_points(O(:,1),2) m_points(O(:,2),2)]');
recx = [ 0 0.2 0.4 0.6 0.8 1  1   1   1   1  1 0.8 0.6 0.4 0.2 0  0   0   0   0  ];
recy = [ 0  0   0   0   0  0 0.2 0.4 0.6 0.8 1  1   1   1   1  1 0.8 0.6 0.4 0.2 ];
rec = [ recx' recy' ]; 

handles.new_border = impoly(gca,rec);
axis 'equal'
%handles.n_b = getPosition(handles.new_border);
guidata(hObject, handles);


% --- Executes on button press in uvGen.
function uvGen_Callback(hObject, eventdata, handles)
% hObject    handle to uvGen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
m_faces = evalin('base','m_faces');
t_points = evalin('base','t_points');

E = sort([m_faces(:,1) m_faces(:,2); m_faces(:,2) m_faces(:,3); m_faces(:,3) m_faces(:,1)]')';
[u,m,n] = unique(E,'rows');
counts = accumarray(n(:), 1);
O = u(counts==1,:);
%%%%%%
OO=O;
bp=zeros(length(OO),2);
bp(1,1) = t_points(OO(1,1),1);
bp(1,2) = t_points(OO(1,1),2);
bp(2,1) = t_points(OO(1,2),1);
bp(2,2) = t_points(OO(1,2),2);
index=3;
first = OO(1,1);
cur = OO(1,2);
OO(1,2) = -1;
while cur~=first
    for i=1:length(OO)
        if cur==OO(i,1)
            bp(index,1) = t_points(OO(i,2),1);
            bp(index,2) = t_points(OO(i,2),2);
            index=index+1;
            cur = OO(i,2);
            OO(i,2)= -1;
            break;
        end
    end
    for i=1:length(OO)
        if cur==OO(i,2)
            bp(index,1) = t_points(OO(i,1),1);
            bp(index,2) = t_points(OO(i,1),2);
            index=index+1;
            cur = OO(i,1);
            OO(i,1)= -1;
            break;
        end
    end
end

handles.plotID = plot(bp(:,1),bp(:,2));
if bp(1,1)==bp(length(bp),1)
    if bp(1,2)==bp(length(bp),2)
        bp = bp(1:length(bp)-1,:);
    end
end

handles.bp = bp;
handles.rotIndices = zeros(1,length(bp));
%handles.new_border = impoly(gca,bp);
axis 'equal'
guidata(hObject, handles);


% --- Executes on button press in flip.
function flip_Callback(hObject, eventdata, handles)
% hObject    handle to flip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.toFlip == 0;
    handles.toFlip = 1;
else
    handles.toFlip = 0;
end
guidata(hObject, handles);


% --- Executes on button press in editB.
function editB_Callback(hObject, eventdata, handles)
% hObject    handle to editB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of editB
set(handles.plotID,'Visible','off')
handles.new_border = impoly(gca,handles.bp);
handles.edit = 1;

guidata(hObject, handles);


% --- Executes on button press in rotationIndicesButton.
function rotationIndicesButton_Callback(hObject, eventdata, handles)
% hObject    handle to rotationIndicesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(handles.plotID,'Visible','off')
hold on
bp = handles.bp;
N = length(bp);
handles.vecOfPlots = zeros(1,N);
for i=1:N
    handles.vecOfPlots(i) = plot(bp(i,1),bp(i,2),'r+','LineWidth',2);
    handles.hcmenu(i)=uicontextmenu;
    item1 = uimenu(handles.hcmenu(i), 'Label', sprintf('%i',handles.rotIndices(i)));
    item2 = uimenu(handles.hcmenu(i), 'Label', '+','callback',{@plus_in,i});
    item3 = uimenu(handles.hcmenu(i), 'Label', '-','callback',{@minus_in,i});
    set(handles.vecOfPlots(i), 'uicontextmenu', handles.hcmenu(i));
end


function plus_in(hObject, handles,i)
myhandles = guidata(gcbo);
myhandles.rotIndices(i)=myhandles.rotIndices(i)+1;

myhandles.vecOfPlots(i)=plot(myhandles.bp(i,1),myhandles.bp(i,2),'r+','LineWidth',2);
    myhandles.hcmenu(i)=uicontextmenu;
    item1 = uimenu(myhandles.hcmenu(i), 'Label', sprintf('%i',myhandles.rotIndices(i)));
    item2 = uimenu(myhandles.hcmenu(i), 'Label', '+','callback',{@plus_in,i});
    item3 = uimenu(myhandles.hcmenu(i), 'Label', '-','callback',{@minus_in,i});
    set(myhandles.vecOfPlots(i), 'uicontextmenu', myhandles.hcmenu(i));
guidata(gcbo,myhandles);

function minus_in(hObject, handles,i)
myhandles = guidata(gcbo);
if myhandles.rotIndices(i)==0
    return
end
myhandles.rotIndices(i)=myhandles.rotIndices(i)-1;

myhandles.vecOfPlots(i)=plot(myhandles.bp(i,1),myhandles.bp(i,2),'r+','LineWidth',2);
    myhandles.hcmenu(i)=uicontextmenu;
    item1 = uimenu(myhandles.hcmenu(i), 'Label', sprintf('%i',myhandles.rotIndices(i)));
    item2 = uimenu(myhandles.hcmenu(i), 'Label', '+','callback',{@plus_in,i});
    item3 = uimenu(myhandles.hcmenu(i), 'Label', '-','callback',{@minus_in,i});
    set(myhandles.vecOfPlots(i), 'uicontextmenu', myhandles.hcmenu(i));
guidata(gcbo,myhandles);
