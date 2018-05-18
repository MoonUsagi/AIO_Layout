function varargout = Layout(varargin)
% LAYOUT MATLAB code for Layout.fig
%      LAYOUT, by itself, creates a new LAYOUT or raises the existing
%      singleton*.
%
%      H = LAYOUT returns the handle to a new LAYOUT or the handle to
%      the existing singleton*.
%
%      LAYOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAYOUT.M with the given input arguments.
%
%      LAYOUT('Property','Value',...) creates a new LAYOUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Layout_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Layout_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Layout

% Last Modified by GUIDE v2.5 08-Mar-2018 16:58:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Layout_OpeningFcn, ...
                   'gui_OutputFcn',  @Layout_OutputFcn, ...
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


% --- Executes just before Layout is made visible.
function Layout_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Layout (see VARARGIN)

% Choose default command line output for Layout
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Layout wait for user response (see UIRESUME)
% uiwait(handles.figure1);

setappdata(handles.figure1,'ListNumber',1);


% --- Outputs from this function are returned to the command line.
function varargout = Layout_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function Input_Callback(hObject, eventdata, handles)
% hObject    handle to Input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%讀檔與選擇檔案類型
[FileName PathName FilterIndex]=uigetfile({'*.jpg';'*.png';'*.bmp';'*.gif';'*.tif';});

if FilterIndex
    str=[PathName FileName];
    OriImg = imread(str);
end
axes(handles.axes1);
imshow(OriImg);

setappdata(handles.figure1,'OriImg',OriImg);
setappdata(handles.figure1,'Input',OriImg);

mImg = getappdata(handles.figure1,'Input');
axes(handles.axes2);
imshow(mImg);

handles.listbox1.String = char(FileName);
% --------------------------------------------------------------------
function Output_Callback(hObject, eventdata, handles)
% hObject    handle to Output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Closebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

if strcmp(get(gcf,'selectiontype'),'open')
    % here you write write code, which you wanna be executed afer double-click
    ListValue = handles.listbox1.Value;
    
    if ListValue == 1
        OriImg = getappdata(handles.figure1,'OriImg');
        axes(handles.axes2);
        imshow(OriImg);
        setappdata(handles.figure1,'Input',OriImg);
    else
        TempMix =  ['Temp',num2str(ListValue)];
        oTi_Image = getappdata(handles.figure1,TempMix);
        axes(handles.axes2);
        imshow(oTi_Image);
        setappdata(handles.figure1,'Input',oTi_Image);
    end
end

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

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ImageProperty_Callback(hObject, eventdata, handles)
% hObject    handle to ImageProperty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Prepocessing_Callback(hObject, eventdata, handles)
% hObject    handle to Prepocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gray_Callback(hObject, eventdata, handles)
% hObject    handle to gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
graymImg = rgb2gray(mImg);
axes(handles.axes3);
imshow(graymImg);
setappdata(handles.figure1,'output',graymImg);
%----G-----
s.a = 1;
setappdata(handles.figure1,'Gnumber',s);
%----------

% --------------------------------------------------------------------
function Binarization_Callback(hObject, eventdata, handles)
% hObject    handle to Binarization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Otsu_Callback(hObject, eventdata, handles)
% hObject    handle to Otsu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
bwmImg = im2bw(mImg);
axes(handles.axes3);
imshow(bwmImg);
setappdata(handles.figure1,'output',bwmImg);
%----G-----
s.a = 2;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Bersen_Callback(hObject, eventdata, handles)
% hObject    handle to Bersen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');

mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
bernsen = mImg;

w = 1;%矩陣大小為2*w+1  
T = 0;%閥值大小  
max = 0;  
min = 0;  
[m,n] = size(bernsen);  
T = zeros(m - 2*w,n - 2*w);  
  
%根據bersen算法計算每個像素點的閥值  
for i = (w + 1):(m - w)  
    for j = (w + 1):(n - w)  
        %預設最大最小值
        max = uint8(bernsen(i,j));  
        min = uint8(bernsen(i,j));
        %掃描周遭 (i+k,i+p) 從-1,0,1
        for k = -w:w  
            for p = -w:w  
                if max < uint8(bernsen(i + k,j + p))  
                    max = uint8(bernsen(i + k,j + p));  
                end  
                if min > uint8(bernsen(i + k,j + p))  
                    min = uint8(bernsen(i + k,j + p));  
                end  
            end  
        end  
        %T 自定義區域閥值
        T(i,j) = 0.5*(max + min);  
    end  
end  
for i = (w + 1):(m - w)  
    for j = (w + 1):(n - w)  
        if bernsen(i,j) > T(i,j)  
            bernsen(i,j) = uint8(255);  
        else  
            bernsen(i,j) = uint8(0);  
        end  
    end  
end  
axes(handles.axes3);
imshow(bernsen);  

setappdata(handles.figure1,'output',bernsen);
% --------------------------------------------------------------------
function I_Bernsen_Callback(hObject, eventdata, handles)
% hObject    handle to I_Bernsen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Input image
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
%% otsu
% level=graythresh(mImg);   %%二值化門檻值
% bwI=im2bw(mImg,level);    %%二值化
%% size of matrix
mImg=double(mImg);
[row,col,dim]=size(mImg);
btI=mImg;
C=zeros(row,col);
G=ones(row,col);
I2=ones(row,col);
T=mImg;

%% set parameters
w=5;          %window size
Tcontrast=15; %neighborhood
sigma=10;     %gaussain
beta=0.9;  %black pixel
apha=0.5;
%% bernsen threshold
Gauss=fspecial('gaussian',[w w],sigma);
Gauss1=imfilter(mImg,Gauss);
Gauss2=imfilter(mImg,Gauss);
for i=w+1:row-w
    for j=w+1:col-w
        for k=-w:w
            sumz=sum(sum(Gauss1(i+k,j+k)));
        end
        Gauss1(i,j)=1/(2*w+1)^2*sumz;
        wI=Gauss1(i-w:i+w,j-w:j+w);                 %% kernel size
        wI2=mImg(i-w:i+w,j-w:j+w);                 %% kernel size
        btI(i,j)=0.5*(max(wI(:))+min(wI(:))); %compute window
        btI2(i,j)=0.5*(max(wI2(:))+min(wI2(:))); %compute window
        T(i,j)=beta*((1-apha)*btI2(i,j)+apha*btI(i,j));
        C(i,j)=max(wI2(:))-min(wI2(:));         %neighborhood
    end
end
  for i=1:row
      for j=1:col
          if C(i,j)< Tcontrast  %% 與相鄰像素非常接近的話
              if T(i,j)>=128
                  BerI_C(i,j)=255;
              else
                  BerI_C(i,j)=0;
              end
          else
              if mImg(i,j)>T(i,j)
                  BerI_C(i,j)=255;
              else
                  BerI_C(i,j)=0;
              end
          end
      end
  end
BerI_C(1:row,1:w)=255;
BerI_C(1:row,col-w:col)=255;
BerI_C(1:w,1:col)=255;
BerI_C(row-w:row,1:col)=255;

axes(handles.axes3);
imshow(BerI_C); 

% axes(handles.axes3);
% H= BerI_C;
% H = histogram(H);
setappdata(handles.figure1,'output',BerI_C);
% --------------------------------------------------------------------
function Niblack_Callback(hObject, eventdata, handles)
% hObject    handle to Niblack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');

mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
 
niblackmImg = mImg;
w = 2;
max =0;  
min =0;  
[m,n]= size(niblackmImg);  
T = zeros(m ,n );  
 
%Niblack掃描點
for i =(w +1):(m - w)  
    for j =(w +1):(n - w)     
        sum =0;
        %掃描周遭 (i+k,i+p) 從-1,0,1
        for k =-w:w  
            for l =-w:w  
                sum = sum + uint32(niblackmImg(i + k,j + l));
            end  
        end  
        average = double(sum)/((2*w+1)*(2*w+1));
        s =0;
        for k =-w:w  
            for l =-w:w  
                s = s +   (uint32(niblackmImg(i + k,j + l))- average)*(uint32(niblackmImg(i + k,j + l))- average);
            end  
        end  
        s= sqrt(double(s)/((2*w+1)*(2*w+1)));
       
        T(i,j)= average +0.2*s;
    end  
end  
for i =  1:m
    for j =1:n
        if niblackmImg(i,j)> T(i,j)  
            niblackmImg(i,j)= uint8(255);  
        else 
            niblackmImg(i,j)= uint8(0);  
        end  
    end  
end
axes(handles.axes3);
imshow(niblackmImg); 

setappdata(handles.figure1,'output',niblackmImg);


% --- Executes on button press in DeleteButton.
function DeleteButton_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ListString = handles.listbox1.String;
ListValue = handles.listbox1.Value;
%暫存
Temp_ListValue = ListValue; 

if ListValue > 1 
    ListValue = ListValue-1;
end

if Temp_ListValue > 1
    ListString(Temp_ListValue,:) = [];
end    
handles.listbox1.Value = ListValue;
handles.listbox1.String = ListString;
TempMix =  ['Temp',num2str(Temp_ListValue)];

setappdata(handles.figure1,TempMix,[])
%h = getappdata(handles.figure1,TempMix)

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(hObject, 'currentpoint');
x = pos(1); y=pos(2);
x1 = x-210;
y1 = y-200;

x2 = x-580;
y2 = y-100;

setappdata(handles.figure1,'axes2_x',x1);
setappdata(handles.figure1,'axes2_y',y1);

setappdata(handles.figure1,'axes3_x',x2);
setappdata(handles.figure1,'axes3_y',y2);

handles.axes2_x.String = ['x :' num2str(x1)];
handles.axes2_y.String = ['y :' num2str(y1)];

handles.axes3_x.String = ['x :' num2str(x2)];
handles.axes3_y.String = ['y :' num2str(y2)];
% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Color_Callback(hObject, eventdata, handles)
% hObject    handle to Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function rgb2hsv_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2hsv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
HsvmImg = rgb2hsv(mImg);
axes(handles.axes3);
imshow(HsvmImg);

axes(handles.axes4);
HistHsvmImg= HsvmImg;
H = histogram(HistHsvmImg);

setappdata(handles.figure1,'output',HsvmImg);

%----G-----
s.a = 6;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function rgb2lab_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2lab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
LabmImg = rgb2lab(mImg);
axes(handles.axes3);
imshow(LabmImg);

axes(handles.axes4);
HistLabmImg= LabmImg;
H = histogram(HistLabmImg);

setappdata(handles.figure1,'output',LabmImg);
%----G-----
s.a = 7;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function rgb2ycbcr_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2ycbcr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
YcbcrmImg = rgb2ycbcr(mImg);
axes(handles.axes3);
imshow(YcbcrmImg);

axes(handles.axes4);
HistYcbcrmImg= YcbcrmImg;
H = histogram(HistYcbcrmImg);

setappdata(handles.figure1,'output',YcbcrmImg);
%----G-----
s.a = 8;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function rgb2xyz_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2xyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
XyzmImg = rgb2xyz(mImg);
axes(handles.axes3);
imshow(XyzmImg);

axes(handles.axes4);
HistXyzmImg= XyzmImg;
H = histogram(HistXyzmImg);

setappdata(handles.figure1,'output',XyzmImg);
%----G-----
s.a = 9;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function rgb2ntsc_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2ntsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
YIQmImg = rgb2ntsc(mImg);
axes(handles.axes3);
imshow(YIQmImg);

axes(handles.axes4);
HistYIQmImg= YIQmImg;
H = histogram(HistYIQmImg);

setappdata(handles.figure1,'output',YIQmImg);
%----G-----
s.a = 10;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function rgb2lin_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
LinmImg = rgb2lin(mImg);
axes(handles.axes3);
imshow(LinmImg);

axes(handles.axes4);
HistLinmImg= LinmImg;
H = histogram(HistLinmImg);

setappdata(handles.figure1,'output',LinmImg);
%----G-----
s.a = 11;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Histogram_Callback(hObject, eventdata, handles)
% hObject    handle to Histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function adjustment_Callback(hObject, eventdata, handles)
% hObject    handle to adjustment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%---------Visible--------
handles.Adj_panel.Visible = 'On';
handles.Sobel_uipanel.Visible = 'Off';
%------------------------
mImg = getappdata(handles.figure1,'Input');
%------------------------
lowin = handles.Lowinslider.Value;
lowout = handles.Lowoutslider.Value;
%------------------------
imadmImg = imadjust(mImg,[lowin lowout],[]);
axes(handles.axes3);
imshow(imadmImg);

axes(handles.axes4);
imadHist= imadmImg;
H = histogram(imadHist);

setappdata(handles.figure1,'output',imadmImg);

%----G-----
s.a = 14;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Equalization_Callback(hObject, eventdata, handles)
% hObject    handle to Equalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
EqmImg = histeq(mImg);
axes(handles.axes3);
imshow(EqmImg);

axes(handles.axes4);
HistEq = EqmImg;
H = histogram(HistEq);

setappdata(handles.figure1,'output',EqmImg);

%----G-----
s.a = 15;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Ahe_Callback(hObject, eventdata, handles)
% hObject    handle to Ahe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
AdaptmImg = adapthisteq(mImg);
axes(handles.axes3);
imshow(AdaptmImg);

axes(handles.axes4);
HistAdapt = AdaptmImg;
H = histogram(HistAdapt);

setappdata(handles.figure1,'output',AdaptmImg);
%----G-----
s.a = 16;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function ImageEnhancement_Callback(hObject, eventdata, handles)
% hObject    handle to ImageEnhancement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Black2HSV_Callback(hObject, eventdata, handles)
% hObject    handle to Black2HSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
I = mImg;

aux=size(I);
M=aux(1);
N=aux(2);
t1=clock;
%%% color difference %%%
maxv=max(I(:,:,1),I(:,:,2));
maxv=max(maxv,I(:,:,3));
minv=min(I(:,:,1),I(:,:,2));
minv=min(minv,I(:,:,3));
ss=double(maxv-minv)./double(maxv); %% saturation
% figure(2);
% imshow(uint8(255*s));
%%% generate weighting functions %%%
W=2.^ss;
Wo=5*W;  %% alpha=4
Wn=2*W;  %% beta=2
%%% generate two exposed images %%%
I=double(I);
Io(:,:,1)=Wo.*I(:,:,1);
Io(:,:,2)=Wo.*I(:,:,2);
Io(:,:,3)=Wo.*I(:,:,3);
In(:,:,1)=Wn.*I(:,:,1);
In(:,:,2)=Wn.*I(:,:,2);
In(:,:,3)=Wn.*I(:,:,3);
Io=uint8(Io);   %% show over-exposd image
% figure(3);
% imshow(Io);
% In=uint8(In);   %% show normal-exposed image
% figure(4);
% imshow(In);
%%% fuse the differt images %%%
I=uint8(I);
v1=double(rgb2gray(I));
W1=exp(-0.5*((v1-128)/128).^2);
v2=double(rgb2gray(Io));
W2=exp(-0.5*((v2-128)/128).^2);
v3=double(rgb2gray(In));
W3=exp(-0.5*((v3-128)/128).^2);
Wt=W1+W2+W3;
I=double(I);
Io=double(Io);
In=double(In);
J(:,:,1)=(W1(:,:).*I(:,:,1)+W2(:,:).*Io(:,:,1)+W3(:,:).*In(:,:,1))./Wt;
J(:,:,2)=(W1(:,:).*I(:,:,2)+W2(:,:).*Io(:,:,2)+W3(:,:).*In(:,:,2))./Wt;
J(:,:,3)=(W1(:,:).*I(:,:,3)+W2(:,:).*Io(:,:,3)+W3(:,:).*In(:,:,3))./Wt;
%%%% show results %%%%
t2=clock;
t2-t1
J=uint8(J);

axes(handles.axes3);
imshow(J); 

axes(handles.axes4);
H= J;
H = histogram(H);

setappdata(handles.figure1,'output',J);
%----G-----
s.a = 12;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Black2Histeq_Callback(hObject, eventdata, handles)
% hObject    handle to Black2Histeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
I = mImg;

I=imresize(I,0.5); %長寬各減半
aux=size(I);
M=aux(1);
N=aux(2);
%t1=clock;
J=rgb2ycbcr(I);
J(:,:,1)=histeq(J(:,:,1)); %% histogram equalization Y
Ie=ycbcr2rgb(J);
%t2=clock;
%t2-t1
v=rgb2gray(Ie);
axis tight %使最大值和最小值範圍一至

axes(handles.axes3);
imshow(Ie); 

axes(handles.axes4);
H= Ie;
H = histogram(H);

setappdata(handles.figure1,'output',Ie);
%----G-----
s.a = 13;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MedianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MedianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
medfilt2mImg = medfilt2(mImg);
axes(handles.axes3);
imshow(medfilt2mImg);

axes(handles.axes4);
Histmedfilt2= medfilt2mImg;
H = histogram(Histmedfilt2);

setappdata(handles.figure1,'output',medfilt2mImg);
%----G-----
s.a = 17;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function EntropyFilter_Callback(hObject, eventdata, handles)
% hObject    handle to EntropyFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
entropymImg = entropyfilt(mImg);
axes(handles.axes3);
imshow(entropymImg,[]);

axes(handles.axes4);
Histentropy= entropymImg;
H = histogram(Histentropy);

setappdata(handles.figure1,'output',entropymImg);
%----G-----
s.a = 18;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function AverageFilter_Callback(hObject, eventdata, handles)
% hObject    handle to AverageFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
h = ones(5,5)/25;
averagemImg = imfilter(mImg,h);
axes(handles.axes3);
imshow(averagemImg);

axes(handles.axes4);
Histaverage= averagemImg;
H = histogram(Histaverage);

setappdata(handles.figure1,'output',averagemImg);
%----G-----
s.a = 19;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function GaussianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to GaussianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
gaussmImg = imgaussfilt(mImg,2);
axes(handles.axes3);
imshow(gaussmImg);

axes(handles.axes4);
Histgauss= gaussmImg;
H = histogram(Histgauss);

setappdata(handles.figure1,'output',gaussmImg);
%----G-----
s.a = 20;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function LaplacianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to LaplacianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
alpha = 0.5;
h = fspecial('laplacian',alpha);
laplacianmImg = imfilter(mImg,h);
axes(handles.axes3);
imshow(laplacianmImg);

axes(handles.axes4);
Histlaplacian= laplacianmImg;
H = histogram(Histlaplacian);

setappdata(handles.figure1,'output',laplacianmImg);
%----G-----
s.a = 21;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function WienerFilter_Callback(hObject, eventdata, handles)
% hObject    handle to WienerFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');

wienermImg = wiener2(mImg,[9 9]);
axes(handles.axes3);
imshow(wienermImg);

axes(handles.axes4);
Histwiener = wienermImg;
H = histogram(Histwiener);

setappdata(handles.figure1,'output',wienermImg);
%----G-----
s.a = 22;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function LowpassFliter_Callback(hObject, eventdata, handles)
% hObject    handle to LowpassFliter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');

b = fir1(30,0.3);
h = ftrans2(b);
Lowpass = imfilter(mImg,h);
axes(handles.axes3);
imshow(Lowpass);

axes(handles.axes4);
Histlowpass = h;
freqz2(Histlowpass);

setappdata(handles.figure1,'output',Lowpass);
%----G-----
s.a = 23;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function HighpassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to HighpassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');

b2 = fir1(30, 0.3, 'high');
h2 = ftrans2(b2);
Highpass = imfilter(mImg, h2, 'symmetric');
axes(handles.axes3);
imshow(Highpass);

axes(handles.axes4);
Histhighpass = h2;
freqz2(Histhighpass);

setappdata(handles.figure1,'output',Highpass);
%----G-----
s.a = 24;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function BandpassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to BandpassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');

b3 = fir1(30, [0.3 0.4]);
h3 = ftrans2(b3);
Bandpass = imfilter(mImg, h3 ,'symmetric');
axes(handles.axes3);
imshow(Bandpass);

axes(handles.axes4);
Histbandpass = h3;
freqz2(Histbandpass);

setappdata(handles.figure1,'output',Bandpass);
%----G-----
s.a = 25;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function ButterworthFilter_Callback(hObject, eventdata, handles)
% hObject    handle to ButterworthFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
fc = 32;
fs = 256;

[b,a] = butter(6,fc/(fs/2));
h4 = ftrans2(b,a);
Butterworth= imfilter(mImg,h4);
axes(handles.axes3);
imshow(Butterworth);

axes(handles.axes4);
Histbutterworth = h4;
freqz2(Histbutterworth);

setappdata(handles.figure1,'output',Butterworth);
%----G-----
s.a = 26;
setappdata(handles.figure1,'Gnumber',s);
%----------

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ocr_Callback(hObject, eventdata, handles)
% hObject    handle to ocr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ocrTrainer;

% --------------------------------------------------------------------
function Hough_Callback(hObject, eventdata, handles)
% hObject    handle to Hough (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Edge_Callback(hObject, eventdata, handles)
% hObject    handle to Edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%-------Visible--------
handles.Sobel_uipanel.Visible = 'On';
handles.Adj_panel.Visible = 'Off';
%----------------------
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
%-----
threshold = 0.05;
%-----
SobelmImg = edge(mImg,'sobel',threshold);
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg);
%----G-----
s.a = 29;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Pewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Pewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
SobelmImg = edge(mImg,'Prewitt');
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg);
%----G-----
s.a = 30;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Roberts_Callback(hObject, eventdata, handles)
% hObject    handle to Roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
SobelmImg = edge(mImg,'Roberts');
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg);
%----G-----
s.a = 31;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Laplacian_Callback(hObject, eventdata, handles)
% hObject    handle to Laplacian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
SobelmImg = edge(mImg,'log');
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg)
%----G-----
s.a = 32;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
SobelmImg = edge(mImg,'Canny');
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg);
%----G-----
s.a = 33;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function Approxcanny_Callback(hObject, eventdata, handles)
% hObject    handle to Approxcanny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
ApproxcannymImg = edge(mImg,'approxcanny');
axes(handles.axes3);
imshow(ApproxcannymImg);

setappdata(handles.figure1,'output',ApproxcannymImg);
%----G-----
s.a = 34;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function FindCircle_Callback(hObject, eventdata, handles)
% hObject    handle to FindCircle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HoughCPanel.Visible = 'On';

mImg = getappdata(handles.figure1,'Input');

min = handles.HC_min_slider.Value;
max = handles.HC_max_slider.Value;

[centers, radii] = imfindcircles(mImg,[min max]);

length(centers)
axes(handles.axes3);
imshow(mImg);

circle = viscircles(centers,radii);
setappdata(handles.figure1,'output',circle);
%----G-----
s.a = 27;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --------------------------------------------------------------------
function FindLine_Callback(hObject, eventdata, handles)
% hObject    handle to FindLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end

rotI = imrotate(mImg,33,'crop');
BW = edge(rotI,'canny');
[H,T,R] = hough(BW);
axes(handles.axes4);
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;

P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
axes(handles.axes3);
plot(x,y,'s','color','white');

lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);

imshow(rotI), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');

setappdata(handles.figure1,'output',rotI);
%setappdata(handles.figure1,'output',circle);
%----G-----
s.a = 28;
setappdata(handles.figure1,'Gnumber',s);
%----------
% --- Executes during object creation, after setting all properties.
function SobelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SobelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function SobelSilder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SobelSilder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function SobelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SobelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SobelEdit as text
%        str2double(get(hObject,'String')) returns contents of SobelEdit as a double

% --- Executes on slider movement.
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
%-----
threshold = str2num(handles.SobelEdit.String)
handles.SobelSilder.Value = threshold;
handles.SobelEdit.String = threshold;
%-----
SobelmImg = edge(mImg,'sobel',threshold);
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg);

function SobelSilder_Callback(hObject, eventdata, handles)
% hObject    handle to SobelSilder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
%-----
threshold = handles.SobelSilder.Value
handles.SobelEdit.Value = threshold;
handles.SobelEdit.String = threshold;
%-----
SobelmImg = edge(mImg,'sobel',threshold);
axes(handles.axes3);
imshow(SobelmImg);

setappdata(handles.figure1,'output',SobelmImg);

% --- Executes during object creation, after setting all properties.
function Lowinslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lowinslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes during object creation, after setting all properties.
function Lowinedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lowinedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function Lowoutslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lowoutslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes during object creation, after setting all properties.
function Lowoutedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lowoutedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on slider movement.
function Lowoutslider_Callback(hObject, eventdata, handles)
% hObject    handle to Lowoutslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mImg = getappdata(handles.figure1,'Input');
%------------------------
lowin = handles.Lowinslider.Value;
lowout = handles.Lowoutslider.Value;

handles.Lowoutedit.Value = lowout;
handles.Lowoutedit.String = lowout;
%------------------------
imadmImg = imadjust(mImg,[lowin lowout],[]);
axes(handles.axes3);
imshow(imadmImg);

axes(handles.axes4);
imadHist= imadmImg;
H = histogram(imadHist);

setappdata(handles.figure1,'output',imadmImg);

function Lowoutedit_Callback(hObject, eventdata, handles)
% hObject    handle to Lowoutedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lowoutedit as text
%        str2double(get(hObject,'String')) returns contents of Lowoutedit as a double
mImg = getappdata(handles.figure1,'Input');
%------------------------
lowin = str2num(handles.Lowinedit.String);
lowout = str2num(handles.Lowoutedit.String);

handles.Lowoutedit.String = lowout;
handles.Lowoutslider.Value = lowout;
%------------------------
imadmImg = imadjust(mImg,[lowin lowout],[]);
axes(handles.axes3);
imshow(imadmImg);

axes(handles.axes4);
imadHist= imadmImg;
H = histogram(imadHist);

setappdata(handles.figure1,'output',imadmImg);

% --- Executes on slider movement.
function Lowinslider_Callback(hObject, eventdata, handles)
% hObject    handle to Lowinslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mImg = getappdata(handles.figure1,'Input');
%------------------------
lowin = handles.Lowinslider.Value;
lowout = handles.Lowoutslider.Value;

handles.Lowinedit.Value = lowin;
handles.Lowinedit.String = lowin;
%------------------------
imadmImg = imadjust(mImg,[lowin lowout],[]);
axes(handles.axes3);
imshow(imadmImg);

axes(handles.axes4);
imadHist= imadmImg;
H = histogram(imadHist);

setappdata(handles.figure1,'output',imadmImg);

function Lowinedit_Callback(hObject, eventdata, handles)
% hObject    handle to Lowinedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Lowinedit as text
%        str2double(get(hObject,'String')) returns contents of Lowinedit as a double
mImg = getappdata(handles.figure1,'Input');
%------------------------
lowin = str2num(handles.Lowinedit.String);
lowout = str2num(handles.Lowoutedit.String);

handles.Lowinedit.String = lowin;
handles.Lowinslider.Value = lowin;
%------------------------
imadmImg = imadjust(mImg,[lowin lowout],[]);
axes(handles.axes3);
imshow(imadmImg);

axes(handles.axes4);
imadHist= imadmImg;
H = histogram(imadHist);

setappdata(handles.figure1,'output',imadmImg)

% --------------------------------------------------------------------
function Closebutton_Callback(hObject, eventdata, handles)
% hObject    handle to Closebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Adj_panel.Visible = 'Off';

% --- Executes on button press in SobelButton.
function SobelButton_Callback(hObject, eventdata, handles)
% hObject    handle to SobelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Sobel_uipanel.Visible = 'Off';

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.HoughCPanel.Visible = 'Off';

% --- Executes during object creation, after setting all properties.
function HC_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HC_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function HC_min_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HC_min_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes during object creation, after setting all properties.
function HC_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HC_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function HC_max_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HC_max_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function HC_min_slider_Callback(hObject, eventdata, handles)
% hObject    handle to HC_min_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mImg = getappdata(handles.figure1,'Input');
%-----------
min = round(handles.HC_min_slider.Value);
max = round(handles.HC_max_slider.Value);

handles.HC_min_edit.Value = min;
handles.HC_min_edit.String = min;
%-----------
[centers, radii] = imfindcircles(mImg,[min max]);

length(centers)
axes(handles.axes3);
imshow(mImg);

circle = viscircles(centers,radii);
setappdata(handles.figure1,'output',circle);

function HC_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to HC_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HC_min_edit as text
%        str2double(get(hObject,'String')) returns contents of HC_min_edit as a double
mImg = getappdata(handles.figure1,'Input');
%-----------
min = handles.HC_min_edit.Value;
max = handles.HC_max_edit.Value;

handles.HC_min_slider.Value = min;
%-----------
[centers, radii] = imfindcircles(mImg,[min max]);

length(centers)
axes(handles.axes3);
imshow(mImg);

circle = viscircles(centers,radii);
setappdata(handles.figure1,'output',circle);

% --- Executes on slider movement.
function HC_max_slider_Callback(hObject, eventdata, handles)
% hObject    handle to HC_max_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mImg = getappdata(handles.figure1,'Input');
%-----------
min = round(handles.HC_min_slider.Value);
max = round(handles.HC_max_slider.Value);

handles.HC_max_edit.Value = max;
handles.HC_max_edit.String = max;
%-----------
[centers, radii] = imfindcircles(mImg,[min max]);

length(centers)
axes(handles.axes3);
imshow(mImg);

circle = viscircles(centers,radii);
setappdata(handles.figure1,'output',circle);

function HC_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to HC_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HC_max_edit as text
%        str2double(get(hObject,'String')) returns contents of HC_max_edit as a double
mImg = getappdata(handles.figure1,'Input');
%-----------
min = handles.HC_min_edit.Value;
max = handles.HC_max_edit.Value;

handles.HC_max_slider.Value = max;
%-----------
[centers, radii] = imfindcircles(mImg,[min max]);

length(centers)
axes(handles.axes3);
imshow(mImg);

circle = viscircles(centers,radii);
setappdata(handles.figure1,'output',circle);
% --------------------------------------------------------------------
function Connected_Callback(hObject, eventdata, handles)
% hObject    handle to Connected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Connected_8_Callback(hObject, eventdata, handles)
% hObject    handle to Connected_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
s = regionprops(mImg,'Centroid');
centroids = cat(1, s.Centroid);
axes(handles.axes3);
%imshow(s);
imshow(mImg);
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off
% --------------------------------------------------------------------
function Connected_c_Callback(hObject, eventdata, handles)
% hObject    handle to Connected_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mImg = getappdata(handles.figure1,'Input');
mysize = size(mImg);
if numel(mysize)>2
    mImg = rgb2gray(mImg);
end
bw = mImg < 100;
imshow(bw)
axes(handles.axes3);
imshow(bw)
stats = regionprops('table',bw,'Centroid',...
    'MajorAxisLength','MinorAxisLength')
centers = stats.Centroid;
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
hold on
viscircles(centers,radii);
hold off

setappdata(handles.figure1,'output',bw);
% --------------------------------------------------------------------
function bwconncomp_8_Callback(hObject, eventdata, handles)
% hObject    handle to bwconncomp_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s.a = 3;
s.b ={25,50,100};
setappdata(handles.figure1,'Gnumber',s);

%OutImage = bwlabel(mImg,8);
%OutImage = bwconncomp(mImg,8);

% --------------------------------------------------------------------
function bwconncomp_4_Callback(hObject, eventdata, handles)
% hObject    handle to bwconncomp_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s.a = 4;
s.b ={22,55,111};
setappdata(handles.figure1,'Gnumber',s);


% --------------------------------------------------------------------
function Morphology_Callback(hObject, eventdata, handles)
% hObject    handle to Morphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Morphology_Dilate_Callback(hObject, eventdata, handles)
% hObject    handle to Morphology_Dilate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----G-----
s.a = 37;
setappdata(handles.figure1,'Gnumber',s);
%----------

%開啟Mophology Shape& Rectange
handles.MorphologyPanel.Visible = 'On';
handles.MorphologyRadioPanel.Visible = 'On';
handles.Panel_Ball.Visible = 'On';
%開啟對應的Radio與Panel
RadioPanel_Dilate = getappdata(handles.figure1,'RadioPanel');
%設置Morphology模式
setappdata(handles.figure1,'Morphology',1);

%判斷RadioPanel_Dilate 是否為空值，不是則開啟對應者
if isempty(RadioPanel_Dilate) == 1
    handles.Ball_M.String = 'Dilate';
    ValueX = getappdata(handles.figure1,'Ballslider_1');
    ValueY = getappdata(handles.figure1,'Ballslider_2');
    if isempty(ValueX) == 1
        ValueX = 5;
    end
    if isempty(ValueY) == 1
        ValueY = 5;
    end
    ValueX = round(ValueX);
    ValueY = round(ValueY);
    se =  offsetstrel('ball',ValueX,ValueY);
    mImg = getappdata(handles.figure1,'Input');
    DilatemImg = imdilate(mImg,se);
    axes(handles.axes3);
    imshow(DilatemImg);
    setappdata(handles.figure1,'output',DilatemImg);
else 
    switch RadioPanel_Dilate
        %Ball的
        case 1
            handles.Ball_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Ballslider_1');
            ValueY = getappdata(handles.figure1,'Ballslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  offsetstrel('ball',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Diamond
        case 2
            handles.Diamond_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Diamondslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('Diamond',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Disk
        case 3
            handles.Disk_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Diskslider_2');
            ValueY = getappdata(handles.figure1,'Diskslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 4;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 4;
            end
            
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
           
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('disk',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Line
        case 4
            handles.Line_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Lineslider_1');
            ValueY = getappdata(handles.figure1,'Lineslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('line',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Octagon
        case 5
            handles.Qctagon_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Octagonslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 3;
            end
            ValueX = 3*(fix(ValueX/3));  
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('octagon',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Cube
        case 6
            handles.Cube_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Cubeslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('cube',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Sphere
        case 7
            handles.Sphere_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Sphereslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('sphere',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Rectange
        case 8 
            handles.Rectange_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Rectangeslider_1');
            ValueY = getappdata(handles.figure1,'Rectangeslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('rectangle',[ValueX ValueY]); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
        %Square
        case 9
            handles.Square_M.String = 'Dilate';
            ValueX = getappdata(handles.figure1,'Squareslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('square',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            DilatemImg = imdilate(mImg,se);
            axes(handles.axes3);
            imshow(DilatemImg);
            setappdata(handles.figure1,'output',DilatemImg);
    end
end

% --------------------------------------------------------------------
function Morphology_Erode_Callback(hObject, eventdata, handles)
% hObject    handle to Morphology_Erode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----G-----
s.a = 38;
setappdata(handles.figure1,'Gnumber',s);
%----------

%開啟Mophology Shape& Rectange
handles.MorphologyPanel.Visible = 'On';
handles.MorphologyRadioPanel.Visible = 'On';
handles.Panel_Ball.Visible = 'On';
%開啟對應的Radio與Panel
RadioPanel_Erode = getappdata(handles.figure1,'RadioPanel');
%設置Morphology模式
setappdata(handles.figure1,'Morphology',2);

%判斷RadioPanel_Erode 是否為空值，不是則開啟對應者
if isempty(RadioPanel_Erode) == 1
    handles.Ball_M.String = 'Erode';
    ValueX = getappdata(handles.figure1,'Ballslider_1');
    ValueY = getappdata(handles.figure1,'Ballslider_2');
    if isempty(ValueX) == 1
       ValueX = 5;
    end
    if isempty(ValueY) == 1
       ValueY = 5;
    end
    ValueX = round(ValueX); 
    ValueY = round(ValueY); 
    se =  offsetstrel('ball',ValueX,ValueY); 
    mImg = getappdata(handles.figure1,'Input');
    ErodemImg = imerode(mImg,se);
    axes(handles.axes3);
    imshow(ErodemImg);
    setappdata(handles.figure1,'output',ErodemImg);
else 
    switch RadioPanel_Erode
        %Ball
        case 1
            handles.Ball_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Ballslider_1');
            ValueY = getappdata(handles.figure1,'Ballslider_2');
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  offsetstrel('ball',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            %Erode
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Diamond
        case 2
            handles.Diamond_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Diamondslider_1');
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            ValueX = round(ValueX); 
            se =  strel('Diamond',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            %Erode
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Disk
        case 3
            handles.Disk_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Diskslider_2');
            ValueY = getappdata(handles.figure1,'Diskslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 4;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 4;
            end
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
           
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('disk',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Line
        case 4
            handles.Line_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Lineslider_1');
            ValueY = getappdata(handles.figure1,'Lineslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('line',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Octagon
        case 5
            handles.Qctagon_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Octagonslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            ValueX = 3*(fix(ValueX/3));  
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('octagon',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Cube
        case 6
            handles.Cube_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Cubeslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('cube',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Sphere
        case 7
            handles.Sphere_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Sphereslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('sphere',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Rectangle
        case 8 
            handles.Rectange_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Rectangeslider_1');
            ValueY = getappdata(handles.figure1,'Rectangeslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('rectangle',[ValueX ValueY]); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
            
        %Square
        case 9
            handles.Square_M.String = 'Erode';
            ValueX = getappdata(handles.figure1,'Squareslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('square',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ErodemImg = imerode(mImg,se);
            axes(handles.axes3);
            imshow(ErodemImg);
            setappdata(handles.figure1,'output',ErodemImg);
    end
end

% --------------------------------------------------------------------
function Morphology_Opening_Callback(hObject, eventdata, handles)
% hObject    handle to Morphology_Opening (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----G-----
s.a = 39;
setappdata(handles.figure1,'Gnumber',s);
%----------

%開啟Mophology Shape& Rectange
handles.MorphologyPanel.Visible = 'On';
handles.MorphologyRadioPanel.Visible = 'On';
handles.Panel_Ball.Visible = 'On';
%開啟對應的Radio與Panel
RadioPanel_Opening = getappdata(handles.figure1,'RadioPanel');
%設置Morphology模式
setappdata(handles.figure1,'Morphology',3);

%判斷RadioPanel_Opening 是否為空值，不是則開啟對應者
if isempty(RadioPanel_Opening) == 1
    handles.Ball_M.String = 'Dilate';
    ValueX = getappdata(handles.figure1,'Ballslider_1');
    ValueY = getappdata(handles.figure1,'Ballslider_2');
    if isempty(ValueX) == 1
       ValueX = 5;
    end
    if isempty(ValueY) == 1
       ValueY = 5;
    end
    ValueX = round(ValueX); 
    ValueY = round(ValueY); 
    se =  offsetstrel('ball',ValueX,ValueY); 
    mImg = getappdata(handles.figure1,'Input');
    OpeningmImg = imopen(mImg,se);
    axes(handles.axes3);
    imshow(OpeningmImg);
    setappdata(handles.figure1,'output',OpeningmImg);
    
else 
    switch RadioPanel_Opening
        %Ball的
        case 1
            handles.Ball_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Ballslider_1');
            ValueY = getappdata(handles.figure1,'Ballslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  offsetstrel('ball',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            %Opening
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Diamond
        case 2
            handles.Diamond_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Diamondslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('Diamond',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Disk
        case 3
            handles.Disk_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Diskslider_2');
            ValueY = getappdata(handles.figure1,'Diskslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 4;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 4;
            end
                        
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
           
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
           
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('disk',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Line
        case 4
            handles.Line_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Lineslider_1');
            ValueY = getappdata(handles.figure1,'Lineslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('line',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Octagon
        case 5
            handles.Qctagon_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Octagonslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            ValueX = 3*(fix(ValueX/3));  
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('octagon',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Cube
        case 6
            handles.Cube_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Cubeslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('cube',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Sphere
        case 7
            handles.Sphere_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Sphereslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('sphere',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Rectangle
        case 8 
            handles.Rectange_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Rectangeslider_1');
            ValueY = getappdata(handles.figure1,'Rectangeslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('rectangle',[ValueX ValueY]); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
            
        %Square
        case 9
            handles.Square_M.String = 'Opening';
            ValueX = getappdata(handles.figure1,'Squareslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('square',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            OpeningmImg = imopen(mImg,se);
            axes(handles.axes3);
            imshow(OpeningmImg);
            setappdata(handles.figure1,'output',OpeningmImg);
    end
end

% --------------------------------------------------------------------
function Morphology_Close_Callback(hObject, eventdata, handles)
% hObject    handle to Morphology_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%----G-----
s.a = 40;
setappdata(handles.figure1,'Gnumber',s);
%----------

%開啟Morphology Shape& Rectange
handles.MorphologyPanel.Visible = 'On';
handles.MorphologyRadioPanel.Visible = 'On';
handles.Panel_Ball.Visible = 'On';
%開啟對應的Radio與Panel
RadioPanel_Close = getappdata(handles.figure1,'RadioPanel');
%設置Morphology模式
setappdata(handles.figure1,'Morphology',4);

%判斷RadioPanel_Dilate 是否為空值，不是則開啟對應者
if isempty(RadioPanel_Close) == 1
    handles.Ball_M.String = 'Close';
    ValueX = getappdata(handles.figure1,'Ballslider_1');
    ValueY = getappdata(handles.figure1,'Ballslider_2');
    if isempty(ValueX) == 1
       ValueX = 5;
    end
    if isempty(ValueY) == 1
       ValueY = 5;
    end
    ValueX = round(ValueX); 
    ValueY = round(ValueY); 
    se =  offsetstrel('ball',ValueX,ValueY); 
    mImg = getappdata(handles.figure1,'Input');
    ClosemImg = imclose(mImg,se);
    axes(handles.axes3);
    imshow(ClosemImg);
    setappdata(handles.figure1,'output',ClosemImg);
    
else 
    switch RadioPanel_Close
        %Ball的
        case 1
            handles.Ball_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Ballslider_1');
            ValueY = getappdata(handles.figure1,'Ballslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  offsetstrel('ball',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            %Close
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Diamond
        case 2
            handles.Diamond_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Diamondslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('Diamond',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Disk
        case 3
            handles.Disk_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Diskslider_2');
            ValueY = getappdata(handles.figure1,'Diskslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 4;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 4;
            end
                        
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
           
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
           
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('disk',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Line
        case 4
            handles.Line_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Lineslider_1');
            ValueY = getappdata(handles.figure1,'Lineslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('line',ValueX,ValueY); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Octagon
        case 5
            handles.Qctagon_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Octagonslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            ValueX = 3*(fix(ValueX/3));  
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('octagon',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Cube
        case 6
            handles.Cube_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Cubeslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('cube',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Sphere
        case 7
            handles.Sphere_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Sphereslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('sphere',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Rectangle
        case 8 
            handles.Rectange_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Rectangeslider_1');
            ValueY = getappdata(handles.figure1,'Rectangeslider_2');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
               ValueX = 5;
            end
            %拉霸Y 空值補預設
            if isempty(ValueY) == 1
               ValueY = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            ValueY = round(ValueY); 
            se =  strel('rectangle',[ValueX ValueY]); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg); 
            setappdata(handles.figure1,'output',ClosemImg);
            
        %Square
        case 9
            handles.Square_M.String = 'Close';
            ValueX = getappdata(handles.figure1,'Squareslider_1');
            %拉霸X 空值補預設
            if isempty(ValueX) == 1
                ValueX = 5;
            end
            %執行運算
            ValueX = round(ValueX); 
            se =  strel('square',ValueX); 
            mImg = getappdata(handles.figure1,'Input');
            ClosemImg = imclose(mImg,se);
            axes(handles.axes3);
            imshow(ClosemImg);
            setappdata(handles.figure1,'output',ClosemImg);
            
     end
end


% --- Executes on button press in rBall.
function rBall_Callback(hObject, eventdata, handles)
% hObject    handle to rBall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rBall

% 關閉其餘Panel
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Ball.Visible = 'On';
% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',1);

MorphologyNumber_r1 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r1
    case 1
          set(handles.Ball_M,'String','Dilate');
    case 2 
          set(handles.Ball_M,'String','Erode');
    case 3 
          set(handles.Ball_M,'String','Opening');
    case 4 
          set(handles.Ball_M,'String','Close');
end

% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',1);


% --- Executes on button press in rDiamond.
function rDiamond_Callback(hObject, eventdata, handles)
% hObject    handle to rDiamond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rDiamond

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Diamond.Visible = 'On';
% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',2);

MorphologyNumber_r2 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r2
    case 1 
          set(handles.Diamond_M,'String','Dilate');
    case 2 
          set(handles.Diamond_M,'String','Erode');
    case 3 
          set(handles.Diamond_M,'String','Opening');
    case 4 
          set(handles.Diamond_M,'String','Close');
end


% --- Executes on button press in rDisk.
function rDisk_Callback(hObject, eventdata, handles)
% hObject    handle to rDisk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rDisk

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';

% 開啟Panel
handles.Panel_Disk.Visible = 'On';
% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',3);

MorphologyNumber_r3 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r3
    case 1 
          set(handles.Disk_M,'String','Dilate');
    case 2 
          set(handles.Disk_M,'String','Erode');
    case 3 
          set(handles.Disk_M,'String','Opening');
    case 4 
          set(handles.Disk_M,'String','Close');
end

% --- Executes on button press in rLine.
function rLine_Callback(hObject, eventdata, handles)
% hObject    handle to rLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rLine

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Line.Visible = 'On';
% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',4);

MorphologyNumber_r4 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r4
    case 1 
          set(handles.Line_M,'String','Dilate');
    case 2 
          set(handles.Line_M,'String','Erode');
    case 3 
          set(handles.Line_M,'String','Opening');
    case 4 
          set(handles.Line_M,'String','Close');
end

% --- Executes on button press in rOctagon.
function rOctagon_Callback(hObject, eventdata, handles)
% hObject    handle to rOctagon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rOctagon

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Octagon.Visible = 'On';
% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',5);

MorphologyNumber_r5 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r5
    case 1 
          set(handles.Qctagon_M,'String','Dilate');
    case 2 
          set(handles.Qctagon_M,'String','Erode');
    case 3 
          set(handles.Qctagon_M,'String','Opening');
    case 4 
          set(handles.Qctagon_M,'String','Close');
end

% --- Executes on button press in rCube.
function rCube_Callback(hObject, eventdata, handles)
% hObject    handle to rCube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rCube

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Cube.Visible = 'On';
% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',6);

MorphologyNumber_r6 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r6
    case 1 
          set(handles.Cube_M,'String','Dilate');
    case 2 
          set(handles.Cube_M,'String','Erode');
    case 3 
          set(handles.Cube_M,'String','Opening');
    case 4 
          set(handles.Cube_M,'String','Close');
end


% --- Executes on button press in rSphere.
function rSphere_Callback(hObject, eventdata, handles)
% hObject    handle to rSphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rSphere

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Sphere.Visible = 'On';

% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',7);

MorphologyNumber_r7 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r7
    case 1 
          set(handles.Sphere_M,'String','Dilate');
    case 2 
          set(handles.Sphere_M,'String','Erode');
    case 3 
          set(handles.Sphere_M,'String','Opening');
    case 4 
          set(handles.Sphere_M,'String','Close');
end

% --- Executes on button press in rRectange.
function rRectange_Callback(hObject, eventdata, handles)
% hObject    handle to rRectange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rRectange

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
handles.Panel_Square.Visible = 'Off';
% 開啟Panel
handles.Panel_Rectange.Visible = 'On';

% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',8);

MorphologyNumber_r8 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r8
    case 1 
          set(handles.Rectange_M,'String','Dilate');
    case 2 
          set(handles.Rectange_M,'String','Erode');
    case 3 
          set(handles.Rectange_M,'String','Opening');
    case 4 
          set(handles.Rectange_M,'String','Close');
end

% --- Executes on button press in rSquare.
function rSquare_Callback(hObject, eventdata, handles)
% hObject    handle to rSquare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rSquare

% 關閉其餘Panel
handles.Panel_Ball.Visible = 'Off';
handles.Panel_Diamond.Visible = 'Off';
handles.Panel_Disk.Visible = 'Off';
handles.Panel_Line.Visible = 'Off';
handles.Panel_Octagon.Visible = 'Off';
handles.Panel_Cube.Visible = 'Off';
handles.Panel_Sphere.Visible = 'Off';
handles.Panel_Rectange.Visible = 'Off';
% 開啟Panel
handles.Panel_Square.Visible = 'On';

% 透過radio 設置開啟對應的Panel
setappdata(handles.figure1,'RadioPanel',9);

MorphologyNumber_r9 = getappdata(handles.figure1,'Morphology');
switch MorphologyNumber_r9
    case 1 
          set(handles.Square_M,'String','Dilate');
    case 2 
          set(handles.Square_M,'String','Erode');
    case 3 
          set(handles.Square_M,'String','Opening');
    case 4 
          set(handles.Square_M,'String','Close');
end


% --- Executes on slider movement.
function Ballslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Ballslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%得到Morphology模式
Morphologynumber =  getappdata(handles.figure1,'Morphology');

%拉霸下的形態學模式
switch Morphologynumber
    %Dilate
    case 1
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
        
    %Erode
    case 2
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %imerode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    %Opening
    case 3
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
     %Close
    case 4
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Ballslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ballslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Ballslider_2_Callback(hObject, eventdata, handles)
% hObject    handle to Ballslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%得到Morphology模式
Morphologynumber =  getappdata(handles.figure1,'Morphology');

%拉霸下的形態學模式
switch Morphologynumber
    %Dilate
    case 1
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
        
    %Erode
    case 2
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
     %Opening
     case 3
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
     %Close
    case 4
        ValueX = handles.Ballslider_1.Value;
        setappdata(handles.figure1,'Ballslider_1',ValueX);
        ValueY = handles.Ballslider_2.Value;
        setappdata(handles.figure1,'Ballslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  offsetstrel('ball',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
end


% --- Executes during object creation, after setting all properties.
function Ballslider_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ballslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Diamondslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Diamondslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Diamond 
    case 1
        ValueX = handles.Diamondslider_1.Value;
        setappdata(handles.figure1,'Diamondslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('diamond',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = handles.Diamondslider_1.Value;
        setappdata(handles.figure1,'Diamondslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('diamond',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
        ValueX = handles.Diamondslider_1.Value;
        setappdata(handles.figure1,'Diamondslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('diamond',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = handles.Diamondslider_1.Value;
        setappdata(handles.figure1,'Diamondslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('diamond',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Diamondslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Diamondslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Diskslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Diskslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Disk 
    case 1
        ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);
        
        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);
        
        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
        ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);
        
        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);
        
        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Diskslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Diskslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Diskslider_2_Callback(hObject, eventdata, handles)
% hObject    handle to Diskslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Diamond 
    case 1
        ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);
        
        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);
        
        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
         ValueX = handles.Diskslider_2.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_2',ValueX);

         ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);
        
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = handles.Diskslider_1.Value;
           if(ValueX >= 0 && ValueX < 4)
              ValueX = 0;
           elseif(ValueX >= 4 && ValueX < 6)
              ValueX = 4;
           elseif(ValueX >= 6 && ValueX < 8)
              ValueX = 6;
           else 
              ValueX = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueX);

        ValueY = handles.Diskslider_1.Value;
           if(ValueY >= 0 && ValueY < 4)
              ValueY = 0;
           elseif(ValueY >= 4 && ValueY < 6)
              ValueY = 4;
           elseif(ValueY >= 6 && ValueY < 8)
                ValueY = 6;
           else 
                ValueY = 8;
           end
        setappdata(handles.figure1,'Diskslider_1',ValueY);

        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('disk',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end


% --- Executes during object creation, after setting all properties.
function Diskslider_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Diskslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Lineslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Lineslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%得到Morphology模式
Morphologynumber =  getappdata(handles.figure1,'Morphology');

%拉霸下的形態學模式
switch Morphologynumber
    %Dilate
    case 1
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
        
    %Erode
    case 2
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %imerode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    %Opening
    case 3
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
     %Close
    case 4
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Lineslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lineslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Lineslider_2_Callback(hObject, eventdata, handles)
% hObject    handle to Lineslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%得到Morphology模式
Morphologynumber =  getappdata(handles.figure1,'Morphology');

%拉霸下的形態學模式
switch Morphologynumber
    %Dilate
    case 1
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
        
    %Erode
    case 2
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %imerode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    %Opening
    case 3
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
     %Close
    case 4
        ValueX = handles.Lineslider_1.Value;
        setappdata(handles.figure1,'Lineslider_1',ValueX);
        ValueY = handles.Lineslider_2.Value;
        setappdata(handles.figure1,'Lineslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Lineslider_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lineslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Octagonslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Octagonslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Diamond 
    case 1
        ValueX = get(handles.Octagonslider_1,'Value');
        ValueX = 3*(fix(ValueX/3));  
        setappdata(handles.figure1,'Octagonslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('octagon',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = get(handles.Octagonslider_1,'Value');
        ValueX = 3*(fix(ValueX/3));  
        setappdata(handles.figure1,'Octagonslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('octagon',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
        ValueX = get(handles.Octagonslider_1,'Value');
        ValueX = 3*(fix(ValueX/3));  
        setappdata(handles.figure1,'Octagonslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('octagon',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = get(handles.Octagonslider_1,'Value');
        ValueX = 3*(fix(ValueX/3));  
        setappdata(handles.figure1,'Octagonslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('octagon',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Octagonslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Octagonslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Cubeslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Cubeslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Diamond 
    case 1
        ValueX = handles.Cubeslider_1.Value;
        setappdata(handles.figure1,'Cubeslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('cube',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = handles.Cubeslider_1.Value;
        setappdata(handles.figure1,'Cubeslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('cube',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
        ValueX = handles.Cubeslider_1.Value;
        setappdata(handles.figure1,'Cubeslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('cube',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = handles.Cubeslider_1.Value;
        setappdata(handles.figure1,'Cubeslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('cube',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Cubeslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cubeslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Sphereslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cubeslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function Sphereslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Panel_Sphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Diamond 
    case 1
        ValueX = handles.Sphereslider_1.Value;
        setappdata(handles.figure1,'Sphereslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('sphere',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = handles.Sphereslider_1.Value;
        setappdata(handles.figure1,'Sphereslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('sphere',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
        ValueX = handles.Sphereslider_1.Value;
        setappdata(handles.figure1,'Sphereslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('sphere',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = handles.Sphereslider_1.Value;
        setappdata(handles.figure1,'Sphereslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('sphere',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Panel_Sphere_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Panel_Sphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Rectangeslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Rectangeslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%得到Morphology模式
Morphologynumber =  getappdata(handles.figure1,'Morphology');

%拉霸下的形態學模式
switch Morphologynumber
    %Dilate
    case 1
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
        
    %Erode
    case 2
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %imerode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
         
    %Opening
    case 3
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        se =  strel('line',ValueX,ValueY);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
     %Close
    case 4
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end

% --- Executes during object creation, after setting all properties.
function Rectangeslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rectangeslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Rectangeslider_2_Callback(hObject, eventdata, handles)
% hObject    handle to Rectangeslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%得到Morphology模式
Morphologynumber =  getappdata(handles.figure1,'Morphology');

%拉霸下的形態學模式
switch Morphologynumber
    %Dilate
    case 1
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
        
    %Erode
    case 2
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %imerode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
          
    %Opening
    case 3
       ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
     %Close
    case 4
        ValueX = handles.Rectangeslider_1.Value;
        setappdata(handles.figure1,'Rectangeslider_1',ValueX);
        ValueY = handles.Rectangeslider_2.Value;
        setappdata(handles.figure1,'Rectangeslider_2',ValueY);
        ValueX = round(ValueX); 
        ValueY = round(ValueY);
        se =  strel('rectangle',[ValueX ValueY]);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
        
end


% --- Executes during object creation, after setting all properties.
function Rectangeslider_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rectangeslider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Squareslider_1_Callback(hObject, eventdata, handles)
% hObject    handle to Squareslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
Morphologynumber =  getappdata(handles.figure1,'Morphology');

switch Morphologynumber
    %Diamond 
    case 1
        ValueX = handles.Squareslider_1.Value;
        setappdata(handles.figure1,'Squareslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('square',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Dilate
        DilatemImg = imdilate(mImg,se);
        axes(handles.axes3);
        imshow(DilatemImg);
        setappdata(handles.figure1,'output',DilatemImg);
    
    case 2
        ValueX = handles.Squareslider_1.Value;
        setappdata(handles.figure1,'Squareslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('square',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Erode
        ErodemImg = imerode(mImg,se);
        axes(handles.axes3);
        imshow(ErodemImg);
        setappdata(handles.figure1,'output',ErodemImg);
        
    case 3
        ValueX = get(handles.Squareslider_1,'Value');
        setappdata(handles.figure1,'Squareslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('square',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Opening
        OpeningmImg = imopen(mImg,se);
        axes(handles.axes3);
        imshow(OpeningmImg);
        setappdata(handles.figure1,'output',OpeningmImg);
        
    case 4
        ValueX = get(handles.Squareslider_1,'Value');
        setappdata(handles.figure1,'Squareslider_1',ValueX);
        ValueX = round(ValueX); 
        se =  strel('square',ValueX);
        mImg = getappdata(handles.figure1,'Input');
        %Close
        ClosemImg = imclose(mImg,se);
        axes(handles.axes3);
        imshow(ClosemImg);
        setappdata(handles.figure1,'output',ClosemImg);
         
end

% --- Executes during object creation, after setting all properties.
function Squareslider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Squareslider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%單次
Gnumber = getappdata(handles.figure1,'Gnumber');
%總
Past_Gnumber = getappdata(handles.figure1,'Past_Gnumber');

All_Gnumber = [Past_Gnumber,Gnumber];

setappdata(handles.figure1,'Past_Gnumber',All_Gnumber);

%-------Visible----------
handles.Sobel_uipanel.Visible = 'Off'; 
handles.Adj_panel.Visible = 'Off';
handles.HoughCPanel.Visible = 'Off';
handles.MorphologyRadioPanel.Visible = 'Off';
handles.Panel_Ball.Visible = 'Off';
handles.MorphologyPanel.Visible = 'Off';

%-------------------------
%取得listbox的String
ListString = get(handles.listbox1,'String');
%取暫存數值ListNumber
ListNumber = getappdata(handles.figure1,'ListNumber');

%ListString{end+1} = ['Export',num2str(ListNumber)];  %cell時使用
ListMix = ['Export',num2str(ListNumber)];
ListOut = strvcat(ListString,ListMix); 
handles.listbox1.String = ListOut;

%ListNumber +1
ListNumber = ListNumber + 1;
%存ListNumber
setappdata(handles.figure1,'ListNumber',ListNumber);

%拿出暫存圖與存入
OutputImage = getappdata(handles.figure1,'output');
TempMix =  ['Temp',num2str(ListNumber)];
setappdata(handles.figure1,TempMix,OutputImage);

setappdata(handles.figure1,'Input',OutputImage);
axes(handles.axes2);
imshow(OutputImage);

% --- Executes on button press in Gcode.
function Gcode_Callback(hObject, eventdata, handles)
% hObject    handle to Gcode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 取出前面疊加的數字與參數
Past_Gnumber = getappdata(handles.figure1,'Past_Gnumber');
% 疊加的總數
G_end = length(Past_Gnumber);
% 空矩陣
Allnumber = [];

for i =1:G_end
    Firstnumber = Past_Gnumber(i).a;
    Allnumber = [Allnumber,Firstnumber]; 
end

strline = [];
str_all = [string("mImg = imread('Input.jpg');")];

for i=1:G_end
    switchin = Allnumber(i);
    switch(switchin)
        case 1
            strline = string("mImg = rgb2gray(mImg);");

        case 2
            strline = string("mImg = im2bw(mImg);");
            
        case 3
            strline =[string("mysize = size(mImg);");
                      string("if numel(mysize)>2");
                      string("   mImg = rgb2gray(mImg);");
                      string("end");
                      string("bernsen = mImg;");
                      string(" ");
                      string("w = 1;");
                      string("T = 0;");
                      string("max = 0;");
                      string("min = 0;");
                      string("[m,n] = size(bernsen);");
                      string("T = zeros(m - 2*w,n - 2*w);");
                      string("for i = (w + 1):(m - w)");
                      string("    for j = (w + 1):(n - w)");
                      string("        max = uint8(bernsen(i,j));");
                      string("        min = uint8(bernsen(i,j));");
                      string("        for k = -w:w ");
                      string("            for p = -w:w");
                      string("                if max < uint8(bernsen(i + k,j + p))");
                      string("                    max = uint8(bernsen(i + k,j + p));");
                      string("                end");
                      string("                if min > uint8(bernsen(i + k,j + p))");
                      string("                    min = uint8(bernsen(i + k,j + p));");
                      string("                end");
                      string("            end");
                      string("        end");
                      string("        T(i,j) = 0.5*(max + min);");
                      string("    end");
                      string("end");
                      string("for i = (w + 1):(m - w)");
                      string("    for j = (w + 1):(n - w)");
                      string("        if bernsen(i,j) > T(i,j)");
                      string("            bernsen(i,j) = uint8(255);");
                      string("        else");
                      string("            bernsen(i,j) = uint8(0);");
                      string("        end");
                      string("    end");
                      string("end");] 
            
        case 4
            strline = [string("mysize = size(mImg);");
                       string("if numel(mysize)>2");
                       string("    mImg = rgb2gray(mImg);");
                       string("end");
                       string(" ");
                       string("mImg=double(mImg);");
                       string("[row,col,dim]=size(mImg);");
                       string("btI=mImg;");
                       string("C=zeros(row,col);");
                       string("G=ones(row,col);");
                       string("I2=ones(row,col);");
                       string("T=mImg;");
                       string(" ");
                       string("w=5;");
                       string("Tcontrast=15;");
                       string("sigma=10;");
                       string("beta=0.9;");
                       string("apha=0.5;");
                       string(" ");
                       string("Gauss=fspecial('gaussian',[w w],sigma);");
                       string("Gauss1=imfilter(mImg,Gauss);");
                       string("Gauss2=imfilter(mImg,Gauss);");
                       string("for i=w+1:row-w");
                       string("    for j=w+1:col-w");
                       string("        for k=-w:w");
                       string("            sumz=sum(sum(Gauss1(i+k,j+k)));");
                       string("        end");
                       string("        Gauss1(i,j)=1/(2*w+1)^2*sumz;");
                       string("        wI=Gauss1(i-w:i+w,j-w:j+w);");
                       string("        wI2=mImg(i-w:i+w,j-w:j+w);");
                       string("         btI(i,j)=0.5*(max(wI(:))+min(wI(:)));");
                       string("        btI2(i,j)=0.5*(max(wI2(:))+min(wI2(:)));");
                       string("        T(i,j)=beta*((1-apha)*btI2(i,j)+apha*btI(i,j));");
                       string("        C(i,j)=max(wI2(:))-min(wI2(:));");
                       string("    end");
                       string("end");
                       string("  for i=1:row");
                       string("      for j=1:col");
                       string("          if C(i,j)< Tcontrast");
                       string("              if T(i,j)>=128");
                       string("                  BerI_C(i,j)=255;");
                       string("              else");
                       string("                  BerI_C(i,j)=0;");
                       string("              end");
                       string("          else");
                       string("              if mImg(i,j)>T(i,j)");
                       string("                  BerI_C(i,j)=255;");
                       string("              else");
                       string("                  BerI_C(i,j)=0;");
                       string("              end");
                       string("          end");
                       string("      end");
                       string("  end");
                       string("BerI_C(1:row,1:w)=255;");
                       string("BerI_C(1:row,col-w:col)=255;");
                       string("BerI_C(1:w,1:col)=255;");
                       string("BerI_C(row-w:row,1:col)=255;");
                       string("mImg = BerI_C;");]
                       
        case 5
            strline = [string("mysize = size(mImg);");
                       string("if numel(mysize)>2");
                       string("    mImg = rgb2gray(mImg);");
                       string("end");
                       string("niblackmImg = mImg;");
                       string("w = 2;");
                       string("max =0;");
                       string("min =0;");
                       string("[m,n]= size(niblackmImg);");
                       string("T = zeros(m ,n );");
                       string(" ");
                       string("for i =(w +1):(m - w)");
                       string("    for j =(w +1):(n - w)");
                       string("        sum =0;");
                       string(" ");
                       string("for k =-w:w");
                       string("            for l =-w:w");
                       string("                sum = sum + uint32(niblackmImg(i + k,j + l));");
                       string("            end");
                       string("        end");
                       string("        average = double(sum)/((2*w+1)*(2*w+1));");
                       string("        s =0;");
                       string("        for k =-w:w");
                       string("            for l =-w:w");
                       string("                s = s +   (uint32(niblackmImg(i + k,j + l))- average)*(uint32(niblackmImg(i + k,j + l))- average);");
                       string("            end");
                       string("        end");
                       string("        s= sqrt(double(s)/((2*w+1)*(2*w+1)));");
                       string(" ");
                       string("        T(i,j)= average +0.2*s;");
                       string("    end");
                       string("end");
                       string("for i =  1:m");
                       string("    for j =1:n");
                       string("        if niblackmImg(i,j)> T(i,j)");
                       string("            niblackmImg(i,j)= uint8(255);");
                       string("        else");
                       string("            niblackmImg(i,j)= uint8(0);");
                       string("        end");
                       string("    end");
                       string("end");]
        case 6
            strline = string("mImg = rgb2hsv(mImg);");
            
        case 7
            strline = string("mImg = rgb2lab(mImg);"); 
            
        case 8
            strline = string("mImg = rgb2ycbcr(mImg);"); 
            
        case 9
            strline = string("mImg = rgb2xyz(mImg);"); 
            
        case 10
            strline = string("mImg = rgb2ntsc(mImg);"); 
            
        case 11
            strline = string("mImg = rgb2lin(mImg);");
            
        case 12
            strline = [string("I = mImg;");
                       string("aux=size(I);");
                       string("M=aux(1);");
                       string("N=aux(2);");
                       string("maxv=max(I(:,:,1),I(:,:,2));");
                       string("maxv=max(maxv,I(:,:,3));");
                       string("minv=min(I(:,:,1),I(:,:,2));");
                       string("minv=min(minv,I(:,:,3));");
                       string("ss=double(maxv-minv)./double(maxv);");
                       string("W=2.^ss;");
                       string("Wo=5*W;");
                       string("Wn=2*W;");
                       string("I=double(I);");
                       string("Io(:,:,1)=Wo.*I(:,:,1);");
                       string("Io(:,:,2)=Wo.*I(:,:,2);");
                       string("Io(:,:,3)=Wo.*I(:,:,3);");
                       string("In(:,:,1)=Wn.*I(:,:,1);");
                       string("Io(:,:,2)=Wn.*I(:,:,2);");
                       string("Io(:,:,3)=Wn.*I(:,:,3);");
                       string("Io=uint8(Io);");
                       string("I=uint8(I);");
                       string("v1=double(rgb2gray(I));");
                       string("W1=exp(-0.5*((v1-128)/128).^2);");
                       string("v2=double(rgb2gray(Io));");
                       string("W2=exp(-0.5*((v2-128)/128).^2);");
                       string("v3=double(rgb2gray(In));");
                       string("W3=exp(-0.5*((v3-128)/128).^2);");
                       string("Wt=W1+W2+W3;");
                       string("I=double(I);");
                       string("Io=double(Io);");
                       string("In=double(In);");
                       string("J(:,:,1)=(W1(:,:).*I(:,:,1)+W2(:,:).*Io(:,:,1)+W3(:,:).*In(:,:,1))./Wt;");
                       string("J(:,:,2)=(W1(:,:).*I(:,:,2)+W2(:,:).*Io(:,:,2)+W3(:,:).*In(:,:,2))./Wt;");
                       string("J(:,:,3)=(W1(:,:).*I(:,:,3)+W2(:,:).*Io(:,:,3)+W3(:,:).*In(:,:,3))./Wt;");
                       string("mImg=uint8(J);")]
        case 13
            strline = [string("mImg = imresize(I,0.5);");
                       string("aux=size(mImg);");
                       string("M=aux(1);");
                       string("M=aux(2);");
                       string("J=rgb2ycbcr(mImg);");
                       string("J(:,:,1)=histeq(J(:,:,1));");
                       string("Ie=ycbcr2rgb(J);");
                       string("mImg=rgb2gray(Ie);");
                       string("axis tight")]
        case 14
            lowin = handles.Lowinslider.Value;
            lowout = handles.Lowoutslider.Value;
            
            strline = [['lowin = ',num2str(lowin)];
                       ['lowout = ',num2str(lowout)];
                       string("mImg = imadjust(mImg,[lowin lowout],[]);")]
           
        case 15
            strline = string("mImg = histeq(mImg);"); 
            
        case 16
            strline = string("mImg = adapthisteq(mImg);"); 
            
        case 17
            strline = string("mImg = medfilt2(mImg);");
            
        case 18
            strline = string("mImg = entropyfilt(mImg);");
            
        case 19
            strline = [string("h = ones(5,5)/25;");
                string("mImg = imfilter(mImg,h);")]
            
        case 20
            strline = string("mImg = imgaussfilt(mImg,2);");
            
        case 21
            strline = [string("alpha = 0.5;");
                string("h = fspecial('laplacian',alpha);");
                string("mImg = imfilter(mImg,h);")]
        case 22
            strline = string("mImg = wiener2(mImg,[9 9]);");
           
        case 23
            strline = [string("b = fir1(30,0.3);");
                string("h = ftrans2(b);");
                 string("mImg = imfilter(mImg,h)");]
            
        case 24
            strline = [string("h2 = ftrans2(b2);");
            string("Highpass = imfilter(mImg, h2, 'symmetric');")]
        
        case 25
            strline = [string("b3 = fir1(30, [0.3 0.4]);");
                string("h3 = ftrans2(b3);");
                string("mImg = imfilter(mImg, h3 ,'symmetri');")]
            
        case 26
            strline = [string("fc = 32;");
                string("fs = 256;");
                string("[b,a] = butter(6,fc/(fs/2));");
                string("h4 = ftrans2(b,a);");
                string("mImg= imfilter(mImg,h4)")]
        
        case 27
            min = handles.HC_min_slider.Value;
            max = handles.HC_max_slider.Value;
            strline = [['min = ',num2str(min)];
                       ['max = ',num2str(max)];
                       string("[centers, radii] = imfindcircles(mImg,[min max]);");
                       string("length(centers);");
                       string("mImg = viscircles(centers,radii);");]
                   
        case 28
            strline = [string("mysize = size(mImg);");
                       string("if numel(mysize)>2");
                       string("    mImg = rgb2gray(mImg);");
                       string("end");
                       string("rotI = imrotate(mImg,33,'crop');");
                       string("BW = edge(rotI,'canny');");
                       string("[H,T,R] = hough(BW);");
                       string("imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');");
                       string("xlabel('\theta'), ylabel('\rho');");
                       string("axis on, axis normal, hold on;");
                       string(" ");
                       string("P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));");
                       string("x = T(P(:,2)); y = R(P(:,1));");
                       string("plot(x,y,'s','color','white');");
                       string(" ");
                       string("lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);");
                       string("imshow(rotI), hold on");
                       string("max_len = 0;");
                       string("for k = 1:length(lines)");
                       string("   xy = [lines(k).point1; lines(k).point2];");
                       string("   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');");
                       string("   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');");
                       string("   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');");
                       string("   len = norm(lines(k).point1 - lines(k).point2);");
                       string("   if ( len > max_len)");
                       string("      max_len = len;");
                       string("      xy_long = xy;");
                       string("   end");
                       string("end");
                       string("plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');");]
                       
        case 29
            threshold = handles.SobelSilder.Value;
      
            strline = [string("mysize = size(mImg)");
                       string("if numel(mysize)>2");
                       string("mImg = rgb2gray(mImg)");
                       string("end");
                       ["threshold = ",num2str(threshold)];
                       string("mImg = edge(mImg,'sobel',threshold)");]
                   
        case 30
            strline = [string("mysize = size(mImg)");
                       string("if numel(mysize)>2");
                       string("mImg = rgb2gray(mImg)");
                       string("end");
                       string("mImg = edge(mImg,'Prewitt')");]
      
        case 31
            strline = [string("mysize = size(mImg)");
                       string("if numel(mysize)>2");
                       string("mImg = rgb2gray(mImg)");
                       string("end");
                       string("mImg = edge(mImg,'Roberts')");]
                        
        case 32
            strline = [string("mysize = size(mImg)");
                       string("if numel(mysize)>2");
                       string("mImg = rgb2gray(mImg)");
                       string("end");
                       string("mImg = edge(mImg,'log')");]
                        
        case 33
            strline = [string("mysize = size(mImg)");
                       string("if numel(mysize)>2");
                       string("mImg = rgb2gray(mImg)");
                       string("end");
                       string("mImg = edge(mImg,'Canny')");]
                        
        case 34
            strline = [string("mysize = size(mImg)");
                       string("if numel(mysize)>2");
                       string("mImg = rgb2gray(mImg)");
                       string("end");
                       string("mImg = edge(mImg,'approxcanny')");]
                        
        case 35
            strline = [string("mysize = size(mImg);");
                       string("if numel(mysize)>2");
                       string("    mImg = rgb2gray(mImg);");
                       string("end");
                       string("s = regionprops(mImg,'Centroid');");
                       string("centroids = cat(1, s.Centroid);");
                       string(" ");
                       string("hold on");
                       string("plot(centroids(:,1),centroids(:,2), 'b*')");
                       string("hold off")]
                       
         case 36
            strline =[string("mysize = size(mImg);");
                      string("if numel(mysize)>2");
                      string("    mImg = rgb2gray(mImg);");
                      string("end");
                      string("bw = mImg < 100;");
                      string("imshow(bw)");
                      string("stats = regionprops('table',bw,'Centroid','MajorAxisLength','MinorAxisLength')");
                      string("centers = stats.Centroid;");
                      string("diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);");
                      string("radii = diameters/2;");
                      string("hold on");
                      string("viscircles(centers,radii);");
                      string("hold off")]
                  
         case 37
            RadioShape = getappdata(handles.figure1,'RadioPanel');
               
                switch RadioShape
                    case 1
                        ValueX = handles.Ballslider_1.Value;
                        ValueY = handles.Ballslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  offsetstrel('ball',ValueX,ValueY);");
                                   string("mImg = imdilate(mImg,se);");]
                               
                    case 2
                        ValueX = handles.Diamondslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('Diamond',ValueX); ");
                                   string("mImg = imdilate(mImg,se);");]
                        
                    case 3
                        ValueX = handles.Diskslider_1.Value;
                        ValueY = handles.Diskslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('disk',ValueX,ValueY);");
                                   string("mImg = imdilate(mImg,se);");]
                               
                    case 4    
                        ValueX = handles.Lineslider_1.Value;
                        ValueY = handles.Lineslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('line',ValueX,ValueY);");
                                   string("mImg = imdilate(mImg,se);");]
                               
                    case 5
                        ValueX = handles.Octagonslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('octagon',ValueX); ");
                                   string("mImg = imdilate(mImg,se);");]
                               
                    case 6 
                        ValueX = handles.Cubeslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('cube',ValueX); ");
                                   string("mImg = imdilate(mImg,se);");]
                               
                    case 7
                        ValueX = handles.Sphereslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('sphere',ValueX); ");
                                   string("mImg = imdilate(mImg,se);");]
                               
                    case 8
                        ValueX = handles.Rectangeslider_1.Value;
                        ValueY = handles.Rectangeslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('rectange',ValueX,ValueY);");
                                   string("mImg = imdilate(mImg,se);");]
                    case 9
                        ValueX = handles.Squareslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('square',ValueX);");
                                   string("mImg = imdilate(mImg,se);");]
                end
                
        case 38
            RadioShape = getappdata(handles.figure1,'RadioPanel');
               
                switch RadioShape
                    case 1
                        ValueX = handles.Ballslider_1.Value;
                        ValueY = handles.Ballslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  offsetstrel('ball',ValueX,ValueY);");
                                   string("mImg = imerode(mImg,se);");]
                               
                    case 2
                        ValueX = handles.Diamondslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('Diamond',ValueX); ");
                                   string("mImg = imerode(mImg,se);");]
                        
                    case 3
                        ValueX = handles.Diskslider_1.Value;
                        ValueY = handles.Diskslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('disk',ValueX,ValueY);");
                                   string("mImg = imerode(mImg,se);");]
                               
                    case 4    
                        ValueX = handles.Lineslider_1.Value;
                        ValueY = handles.Lineslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('line',ValueX,ValueY);");
                                   string("mImg = imerode(mImg,se);");]
                               
                    case 5
                        ValueX = handles.Octagonslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('octagon',ValueX); ");
                                   string("mImg = imerode(mImg,se);");]
                               
                    case 6 
                        ValueX = handles.Cubeslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('cube',ValueX); ");
                                   string("mImg = imerode(mImg,se);");]
                               
                    case 7
                        ValueX = handles.Sphereslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('sphere',ValueX); ");
                                   string("mImg = imerode(mImg,se);");]
                               
                    case 8
                        ValueX = handles.Rectangeslider_1.Value;
                        ValueY = handles.Rectangeslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('rectange',ValueX,ValueY);");
                                   string("mImg = imerode(mImg,se);");]
                    case 9
                        ValueX = handles.Squareslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('square',ValueX); ");
                                   string("mImg = imerode(mImg,se);");]
                end
                
        case 39
            RadioShape = getappdata(handles.figure1,'RadioPanel');
               
                switch RadioShape
                    case 1
                        ValueX = handles.Ballslider_1.Value;
                        ValueY = handles.Ballslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  offsetstrel('ball',ValueX,ValueY);");
                                   string("mImg = imopen(mImg,se);");]
                               
                    case 2
                        ValueX = handles.Diamondslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('Diamond',ValueX); ");
                                   string("mImg = imopen(mImg,se);");]
                        
                    case 3
                        ValueX = handles.Diskslider_1.Value;
                        ValueY = handles.Diskslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('disk',ValueX,ValueY);");
                                   string("mImg = imopen(mImg,se);");]
                               
                    case 4    
                        ValueX = handles.Lineslider_1.Value;
                        ValueY = handles.Lineslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('line',ValueX,ValueY);");
                                   string("mImg = imopen(mImg,se);");]
                               
                    case 5
                        ValueX = handles.Octagonslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('octagon',ValueX); ");
                                   string("mImg = imopen(mImg,se);");]
                               
                    case 6 
                        ValueX = handles.Cubeslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('cube',ValueX); ");
                                   string("mImg = imopen(mImg,se);");]
                               
                    case 7
                        ValueX = handles.Sphereslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('sphere',ValueX); ");
                                   string("mImg = imopen(mImg,se);");]
                               
                    case 8
                        ValueX = handles.Rectangeslider_1.Value;
                        ValueY = handles.Rectangeslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('rectange',ValueX,ValueY);");
                                   string("mImg = imopen(mImg,se);");]
                    case 9
                        ValueX = handles.Squareslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('square',ValueX); ");
                                   string("mImg = imopen(mImg,se);");]
                end
                
        case 40
            RadioShape = getappdata(handles.figure1,'RadioPanel');
               
                switch RadioShape
                    case 1
                        ValueX = handles.Ballslider_1.Value;
                        ValueY = handles.Ballslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  offsetstrel('ball',ValueX,ValueY);");
                                   string("mImg = imclose(mImg,se);");]
                               
                    case 2
                        ValueX = handles.Diamondslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('Diamond',ValueX); ");
                                   string("mImg = imclose(mImg,se);");]
                        
                    case 3
                        ValueX = handle.Diskslider_1.Value;
                        ValueY = handles.Diskslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('disk',ValueX,ValueY);");
                                   string("mImg = imclose(mImg,se);");]
                               
                    case 4    
                        ValueX = handles.Lineslider_1.Value;
                        ValueY = handles.Lineslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('line',ValueX,ValueY);");
                                   string("mImg = imclose(mImg,se);");]
                               
                    case 5
                        ValueX = handles.Octagonslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('octagon',ValueX); ");
                                   string("mImg = imclose(mImg,se);");]
                               
                    case 6 
                        ValueX = handles.Cubeslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('cube',ValueX); ");
                                   string("mImg = imclose(mImg,se);");]
                               
                    case 7
                        ValueX = handles.Sphereslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('sphere',ValueX); ");
                                   string("mImg = imclose(mImg,se);");]
                               
                    case 8
                        ValueX = handles.Rectangeslider_1.Value;
                        ValueY = handles.Rectangeslider_2.Value;
                        
                        strline = [['ValueX = ',num2str(ValueX)];
                                   ['ValueY = ',num2str(ValueY)];
                                   string("ValueX = round(ValueX);");
                                   string("ValueY = round(ValueY);");
                                   string("se =  strel('rectange',ValueX,ValueY);");
                                   string("mImg = imclose(mImg,se);");]
                    case 9
                        ValueX = handle.Squareslider_1.Value;
                        strline = [['ValueX = ',num2str(ValueX)];
                                   string("ValueX = round(ValueX);");
                                   string("se =  strel('square',ValueX);");
                                   string("mImg = imclose(mImg,se);");]
                end
                
        case 41
            disp('2--222')
            strline = [string("bwmImg = im2bw(mImg);");
                string("bwmImg = im2bw(mImg);kkk2;")];
            %eval(['strline',num2str(i), '=strline']);
    end
    str_all = [str_all;strline];
    disp(str_all);
end
str_end = [string("figure;");
           string("imshow(mImg);")]
       
str_all = [str_all;str_end];

fileID = fopen('Gcode.m','w');

fprintf(fileID,'%s\n',str_all);

fclose(fileID);



function Ball_Edit1_Callback(hObject, eventdata, handles)
% hObject    handle to Ball_Edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ball_Edit1 as text
%        str2double(get(hObject,'String')) returns contents of Ball_Edit1 as a double


% --- Executes during object creation, after setting all properties.
function Ball_Edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ball_Edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ball_Edit2_Callback(hObject, eventdata, handles)
% hObject    handle to Ball_Edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ball_Edit2 as text
%        str2double(get(hObject,'String')) returns contents of Ball_Edit2 as a double


% --- Executes during object creation, after setting all properties.
function Ball_Edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ball_Edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Calibrator_Callback(hObject, eventdata, handles)
% hObject    handle to Calibrator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Camera_Callback(hObject, eventdata, handles)
% hObject    handle to Camera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cameraCalibrator;

% --------------------------------------------------------------------
function StereoCamera_Callback(hObject, eventdata, handles)
% hObject    handle to StereoCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stereoCameraCalibrator;


% --------------------------------------------------------------------
function APP_Callback(hObject, eventdata, handles)
% hObject    handle to APP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Viewer_Callback(hObject, eventdata, handles)
% hObject    handle to Viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function DIcom_Viewer_Callback(hObject, eventdata, handles)
% hObject    handle to DIcom_Viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dicomBrowser();

% --------------------------------------------------------------------
function Image_Browser_Callback(hObject, eventdata, handles)
% hObject    handle to Image_Browser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageBrowser();

% --------------------------------------------------------------------
function Image_Viewer_Callback(hObject, eventdata, handles)
% hObject    handle to Image_Viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imtool();

% --------------------------------------------------------------------
function Map_Viewer_Callback(hObject, eventdata, handles)
% hObject    handle to Map_Viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mapview();

% --------------------------------------------------------------------
function Video_Viewer_Callback(hObject, eventdata, handles)
% hObject    handle to Video_Viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
implay();

% --------------------------------------------------------------------
function Volume_Viewer_Callback(hObject, eventdata, handles)
% hObject    handle to Volume_Viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
volumeViewer();

% --------------------------------------------------------------------
function GroundTruth_Labeler_Callback(hObject, eventdata, handles)
% hObject    handle to GroundTruth_Labeler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
groundTruthLabeler();

% --------------------------------------------------------------------
function Untitled_19_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imagelabeler();

% --------------------------------------------------------------------
function Color_thresholder_Callback(hObject, eventdata, handles)
% hObject    handle to Color_thresholder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colorThresholder();

% --------------------------------------------------------------------
function Image_Acquisition_Callback(hObject, eventdata, handles)
% hObject    handle to Image_Acquisition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imaqtool();

% --------------------------------------------------------------------
function Image_RegionAnalyzer_Callback(hObject, eventdata, handles)
% hObject    handle to Image_RegionAnalyzer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageRegionAnalyzer();

% --------------------------------------------------------------------
function OCRTrainer_Callback(hObject, eventdata, handles)
% hObject    handle to OCRTrainer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ocrTrainer();

% --------------------------------------------------------------------
function RegistrationEstimator_Callback(hObject, eventdata, handles)
% hObject    handle to RegistrationEstimator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
registrationEstimator();


% --------------------------------------------------------------------
function Image_BatchProcessor_Callback(hObject, eventdata, handles)
% hObject    handle to Image_BatchProcessor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageBatchProcessor();

% --------------------------------------------------------------------
function Image_Segmenter_Callback(hObject, eventdata, handles)
% hObject    handle to Image_Segmenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageSegmenter();
