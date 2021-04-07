function varargout = edgeDetectionG(varargin)
% EDGEDETECTIONG MATLAB code for edgeDetectionG.fig
%      EDGEDETECTIONG, by itself, creates a new EDGEDETECTIONG or raises the existing
%      singleton*.
%
%      H = EDGEDETECTIONG returns the handle to a new EDGEDETECTIONG or the handle to
%      the existing singleton*.
%
%      EDGEDETECTIONG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDGEDETECTIONG.M with the given input arguments.
%
%      EDGEDETECTIONG('Property','Value',...) creates a new EDGEDETECTIONG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before edgeDetectionG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to edgeDetectionG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help edgeDetectionG

% Last Modified by GUIDE v2.5 16-Mar-2021 19:37:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @edgeDetectionG_OpeningFcn, ...
                   'gui_OutputFcn',  @edgeDetectionG_OutputFcn, ...
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


% --- Executes just before edgeDetectionG is made visible.
function edgeDetectionG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to edgeDetectionG (see VARARGIN)

% Choose default command line output for edgeDetectionG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes edgeDetectionG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = edgeDetectionG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadImage.
function loadImage_Callback(hObject, eventdata, handles)
% hObject    handle to loadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%img_R = imread('images.jpeg');
% use axes(handles.axes1) to display in axes; 
%imshow(img_R);


[filename, pathname] = uigetfile('*.*', 'Pick a MATLAB code file');
if isequal(filename,0) || isequal(pathname,0)
   disp('You pressed cancel')
else
   filename = strcat(pathname, filename);
   img_R = imread(filename);
   axes(handles.axes1);
   imshow(img_R);
end
handles.img_R = img_R;

%Updating handles Structure
guidata(hObject, handles);

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
img_R = handles.img_R;
value = get(hObject,'Value');
switch value
    case 1
        img = rgb2gray(img_R);
        p_x=[-1 0 1;-1 0 1;-1 0 1];
        imgx=filter2(p_x,img);
        p_y=p_x';
        imgy=filter2(p_y,img);
        edge = sqrt(imgx.^2+imgy.^2);
        edge_thresholding = imbinarize(edge/255,0.8);
        axes(handles.axes2);
        imshow(edge_thresholding);
    case 2
        Nx1=10;Sigmax1=2;Nx2=10;Sigmax2=4;Theta1=pi/2;
        Ny1=10;Sigmay1=2;Ny2=10;Sigmay2=4;Theta2=pi/2;
        alfa=0.9;
        img = rgb2gray(img_R);
        filterx=d2dgauss(Nx1,Sigmax1,Nx2,Sigmax2,Theta1);
        Ix= conv2(img,filterx,'same');
        filtery=d2dgauss(Ny1,Sigmay1,Ny2,Sigmay2,Theta2);
        Iy=conv2(img,filtery,'same'); 
        NVI=sqrt(Ix.*Ix+Iy.*Iy);
        I_max=max(max(NVI));
        I_min=min(min(NVI));
        level=alfa*(I_max-I_min)+I_min;
        
        Ibw=max(NVI,level.*ones(size(NVI)));
        axes(handles.axes2);
        %Ibw = imadd(Ibw, 1000);
        %Ibw = mat2gray(Ibw);
        imshow(Ibw, []);
        %colormap('gray');
end
function h = d2dgauss(n1,sigma1,n2,sigma2,theta)
r=[cos(theta) -sin(theta);
   sin(theta)  cos(theta)];
for i = 1 : n2 
    for j = 1 : n1
        u = r * [j-(n1+1)/2 i-(n2+1)/2]';
        h = zeros(size(u));
        h(i,j) = gauss(u(1),sigma1)*dgauss(u(2),sigma2);
    end
end
h = h / sqrt(sum(sum(abs(h).*abs(h))));

% Function "gauss.m":
function y = gauss(x,std)
y = exp(-x^2/(2*std^2)) / (std*sqrt(2*pi));

% Function "dgauss.m"(first order derivative of gauss function):
function y = dgauss(x,std)
y = -x * gauss(x,std) / std^2;

 
% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
