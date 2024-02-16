function varargout = SugarCaneLeafDisease(varargin)
% SUGARCANELEAFDISEASE MATLAB code for SugarCaneLeafDisease.fig
%      SUGARCANELEAFDISEASE, by itself, creates a new SUGARCANELEAFDISEASE or raises the existing
%      singleton*.
%
%      H = SUGARCANELEAFDISEASE returns the handle to a new SUGARCANELEAFDISEASE or the handle to
%      the existing singleton*.
%
%      SUGARCANELEAFDISEASE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUGARCANELEAFDISEASE.M with the given input arguments.
%
%      SUGARCANELEAFDISEASE('Property','Value',...) creates a new SUGARCANELEAFDISEASE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SugarCaneLeafDisease_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SugarCaneLeafDisease_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SugarCaneLeafDisease

% Last Modified by GUIDE v2.5 26-May-2017 00:09:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SugarCaneLeafDisease_OpeningFcn, ...
                   'gui_OutputFcn',  @SugarCaneLeafDisease_OutputFcn, ...
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


% --- Executes just before SugarCaneLeafDisease is made visible.
function SugarCaneLeafDisease_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SugarCaneLeafDisease (see VARARGIN)

% Choose default command line output for SugarCaneLeafDisease
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SugarCaneLeafDisease wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SugarCaneLeafDisease_OutputFcn(hObject, eventdata, handles) 
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
clc;
[filename, pathname] = uigetfile({'*.*';'*.bmp';'*.jpg';'*.gif'}, 'Pick a Leaf Image File');
I = imread([pathname,filename]);
I = imresize(I,[256,256]);

axes(handles.axes1);
imshow(I);title('ORIGINAL IMAGE');
handles.ImgData1 = I;
guidata(hObject,handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.ImgData1;
I = rgb2gray(I);
%figure,imshow(I);

I = imresize(I,[256,256]);
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
axes(handles.axes2);
imshow(gradmag,[]), title('Gradient magnitude (gradmag)')


handles.ImgData2 = gradmag;
guidata(hObject,handles);
%figure
%imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
%gradmag=~gradmag;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.ImgData1;
I=imresize(I,[256 256]);

%lab_he=rgb2hsv(I);
cform = makecform('srgb2lab');
lab_he = applycform(I,cform);
axes(handles.axes3);
imshow(lab_he),title('L*a*b');
handles.ImgData3 = lab_he;
guidata(hObject,handles);





% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imge=handles.ImgData2;
imge=imresize(imge,[256 256]);
I=imge;
[ T ] = otsurec(I, 8);

ttotal=8;
[mt]=size(T);

I1=I;
[N_St] = Normalizefunction(I1);
figure,
for i=1:mt
    if i==ttotal
        t2=1;
        %t2
    else
        t2=T(i+1);
        %t2
    end
    t1=T(i);
    %t1
    [ Ib ] = binimage(N_St,t1,t2);
    %Store the resultent bimary image in the following folder, before tht
    %you have to create that folder....
    %dest='D:\An efficient algorithm for fractal analysis of textures\image\';
    fs1=strcat(int2str(i),'.bmp');
    imwrite(Ib, fs1);
    subplot(2,4,i),imshow(Ib),title('Otsu Threshold Image');
end
x = inputdlg('Enter Otsu threshold image number {4 or 6 or 8}:');
k = str2num(x{:});

if k==4
    x=4;
    inimg=strcat(int2str(x),'.bmp');
img=imread(inimg);
img1=~img;
axes(handles.axes5);
imshow(img1),title('Otsu Thresholding [Detected] Image');
handles.ImgData4 = img1;
guidata(hObject,handles);
end
if k==6
    x=6;
    inimg=strcat(int2str(x),'.bmp');
img=imread(inimg);
img2=~img;
axes(handles.axes5);
imshow(img2),title('Otsu Thresholding [Detected] Image');
handles.ImgData4 = img2;
guidata(hObject,handles);
end
if k==8
    x=8;
    inimg=strcat(int2str(x),'.bmp');
img3=imread(inimg);
axes(handles.axes5);
imshow(img3),title('Otsu Thresholding [Detected] Image');
handles.ImgData4 = img3;
guidata(hObject,handles);
end


x=6;


x=8;


axes(handles.axes5);
imshow(img3),title('Otsu Thresholding [Detected] Image');
handles.ImgData4 = img3;
guidata(hObject,handles);




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lab_he=handles.ImgData3;
lab_he=imresize(lab_he,[256 256]);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

I=handles.ImgData1;
I=imresize(I,[256 256]);

nColors = 3;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);

pixel_labels = reshape(cluster_idx,nrows,ncols);
figure,imshow(pixel_labels,[]), title('image labeled by cluster index');

segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = I;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
    fs1=strcat(int2str(k),'.bmp');
    imwrite(segmented_images{k}, fs1);
end
axes(handles.axes4);
imshow(segmented_images{1}), title('objects in cluster 1');
handles.ImgData51 = segmented_images{1};
guidata(hObject,handles);
axes(handles.axes7);
imshow(segmented_images{2}), title('objects in cluster 2');
handles.ImgData52 = segmented_images{2};
guidata(hObject,handles);
axes(handles.axes8);
imshow(segmented_images{3}), title('objects in cluster 3');
handles.ImgData53 = segmented_images{3};
guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img3=handles.ImgData4;
img=imresize(img3,[256 256]);

s=strel('disk',2,0);
I16=imerode(img,s);
I17=(img-I16);
%I17=imclose(I17,s);
%I18 = imresize(I17,[300,400]);
%figure,imshow(I17);title('BOUNDARY EXTRACTED');

AVG_NFT=I17
I=handles.ImgData1;
A=imresize(I,[256 256]);

bin = im2bw(A,0.5);
[m,n]=size(AVG_NFT);
MFC1=0;
for i=1:m
    for j=1:n
        
        MFC1=MFC1+double(AVG_NFT(i,j));
        
    end
end
MFC=(MFC1/(m*n));

count=0;
for i=1:m
    for j=1:n
        
        if AVG_NFT(i,j) >= MFC
            count=count+1;
        end
    end
end

%M = FCM_cluster(A);
%M =~im2bw(M,.5);
%figure,imshow(M);title('clustered image');
%bin=M;


%bin=Kmeansclustering(A);


SFC=MFC1;
NHF=count;
T=(SFC/((m*n)+NHF));

T=uint8(T);
for i=1:m
    for j=1:n
        
        if (AVG_NFT(i,j)>T && bin(i,j)==1)
            bin(i,j)=1;
            res(i,j)=1;
        else
            bin(i,j)=0;
            res(i,j)=0;
        end
        
    end
end

B=res;
bin=B;



I13=bin;
I14 = imresize(I13,[256,256]);
%figure,
axes(handles.axes6);
imshow(I14);title(' THRESHOLDED IMAGE');
handles.ImgData6 = I14;
guidata(hObject,handles);




% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = inputdlg('Enter Cluster number {1 or 2 or 3}:');
k = str2num(x{:}); 
if k==1
Cluster_selected=handles.ImgData51;
end
if k==2
        Cluster_selected=handles.ImgData52;
    end
   if k==3
        Cluster_selected=handles.ImgData53;
   end

        g1=rgb2gray(Cluster_selected);
         %figure,imshow(g1),title('ROI [Extracted Image]');
    
         
            A=g1;
            AVG_NFT=g1;
            bin = im2bw(A,0.5);
           % figure,imshow(bin);title('Bin T=0.5');
            [m,n]=size(AVG_NFT);
            MFC1=0;
            for i=1:m
                for j=1:n
                    MFC1=MFC1+double(AVG_NFT(i,j));
                end
            end
            MFC=(MFC1/(m*n));

            count=0;
            for i=1:m
                for j=1:n

                    if AVG_NFT(i,j) >= MFC
                        count=count+1;
                    end
                end
            end

            SFC=MFC1;
            NHF=count;
            T=(SFC/((m*n)+NHF));

            T=uint8(T);
            for i=1:m
                for j=1:n

                    if (AVG_NFT(i,j)>T && bin(i,j)==1)
                        bin(i,j)=1;
                        res(i,j)=1;
                    else
                        bin(i,j)=0;
                        res(i,j)=0;
                    end

                end
            end

            B=res;
   
   
g1=rgb2gray(Cluster_selected);
e1=edge(g1,'sobel');
axes(handles.axes9);
imshow(e1),title('ROI [Extracted Image]');
handles.ImgData7 = e1;
handles.ImgDataB = B;
guidata(hObject,handles);



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I19=handles.ImgData6;
glcms = graycomatrix(I19);
% Derive Statistics from GLCM
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast = stats.Contrast;
Correlation = stats.Correlation;
Energy = stats.Energy;
Homogeneity = stats.Homogeneity;
Mean = mean2(I19);
Standard_Deviation = std2(I19);
Entropy = entropy(I19);
RMS = mean2(rms(I19));
Variance = mean2(var(double(I19)));
a = sum(double(I19(:)));
Smoothness = 1-(1/(1+a));
Kurtosis = kurtosis(double(I19(:)));
Skewness = skewness(double(I19(:)));
% Inverse Difference Movement
m = size(I19,1);
n = size(I19,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = I19(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM = double(in_diff);
feat_disease = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, RMS, Variance, Smoothness, Kurtosis, Skewness, IDM];
set(handles.edit1,'string',Contrast);
set(handles.edit2,'string',Correlation);
set(handles.edit3,'string',Energy);
set(handles.edit4,'string',Homogeneity);
set(handles.edit5,'string',Mean);
set(handles.edit6,'string',Standard_Deviation);
set(handles.edit7,'string',Entropy);
set(handles.edit8,'string',RMS);
set(handles.edit9,'string',Variance);
set(handles.edit10,'string',Smoothness);
set(handles.edit11,'string',Kurtosis);
set(handles.edit12,'string',Skewness);
set(handles.edit13,'string',IDM);
guidata(hObject,handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I19=handles.ImgData7;
glcms = graycomatrix(I19);
% Derive Statistics from GLCM
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast = stats.Contrast;
Correlation = stats.Correlation;
Energy = stats.Energy;
Homogeneity = stats.Homogeneity;
Mean = mean2(I19);
Standard_Deviation = std2(I19);
Entropy = entropy(I19);
%rms(I19);
%RMS = mean2(rms(I19));
Variance = mean2(var(double(I19)));
a = sum(double(I19(:)));
Smoothness = 1-(1/(1+a));
Kurtosis = kurtosis(double(I19(:)));
Skewness = skewness(double(I19(:)));
% Inverse Difference Movement
m = size(I19,1);
n = size(I19,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = I19(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM = double(in_diff);

I19=handles.ImgDataB;
glcms = graycomatrix(I19);
% Derive Statistics from GLCM
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast1 = stats.Contrast;
Correlation1 = stats.Correlation;
Energy1 = stats.Energy;
Homogeneity1 = stats.Homogeneity;
Mean1 = mean2(I19);
Standard_Deviation1 = std2(I19);
Entropy1 = entropy(I19);
%rms(I19);
%RMS = mean2(rms(I19));
Variance1 = mean2(var(double(I19)));
a = sum(double(I19(:)));
Smoothness1 = 1-(1/(1+a));
Kurtosis1 = kurtosis(double(I19(:)));
Skewness1 = skewness(double(I19(:)));
% Inverse Difference Movement
m = size(I19,1);
n = size(I19,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = I19(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM1 = double(in_diff);
feat_color_disease = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, Variance, Smoothness, Kurtosis, Skewness, IDM,Contrast1,Correlation1,Energy1,Homogeneity1, Mean1, Standard_Deviation1, Entropy1, Variance1, Smoothness1, Kurtosis1, Skewness1, IDM1];

set(handles.edit14,'string',Contrast);
set(handles.edit15,'string',Correlation);
set(handles.edit16,'string',Energy);
set(handles.edit17,'string',Homogeneity);
set(handles.edit18,'string',Mean);
set(handles.edit19,'string',Standard_Deviation);
set(handles.edit20,'string',Entropy);

set(handles.edit22,'string',Variance);
set(handles.edit23,'string',Smoothness);
set(handles.edit24,'string',Kurtosis);
set(handles.edit25,'string',Skewness);
set(handles.edit26,'string',IDM);
guidata(hObject,handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I19=handles.ImgData7;
glcms = graycomatrix(I19);
% Derive Statistics from GLCM
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast = stats.Contrast;
Correlation = stats.Correlation;
Energy = stats.Energy;
Homogeneity = stats.Homogeneity;
Mean = mean2(I19);
Standard_Deviation = std2(I19);
Entropy = entropy(I19);
%rms(I19);
%RMS = mean2(rms(I19));
Variance = mean2(var(double(I19)));
a = sum(double(I19(:)));
Smoothness = 1-(1/(1+a));
Kurtosis = kurtosis(double(I19(:)));
Skewness = skewness(double(I19(:)));
% Inverse Difference Movement
m = size(I19,1);
n = size(I19,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = I19(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM = double(in_diff);

I19=handles.ImgDataB;
glcms = graycomatrix(I19);
% Derive Statistics from GLCM
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast1 = stats.Contrast;
Correlation1 = stats.Correlation;
Energy1 = stats.Energy;
Homogeneity1 = stats.Homogeneity;
Mean1 = mean2(I19);
Standard_Deviation1 = std2(I19);
Entropy1 = entropy(I19);
%rms(I19);
%RMS = mean2(rms(I19));
Variance1 = mean2(var(double(I19)));
a = sum(double(I19(:)));
Smoothness1 = 1-(1/(1+a));
Kurtosis1 = kurtosis(double(I19(:)));
Skewness1 = skewness(double(I19(:)));
% Inverse Difference Movement
m = size(I19,1);
n = size(I19,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = I19(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM1 = double(in_diff);
feat_color_disease = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, Variance, Smoothness, Kurtosis, Skewness, IDM,Contrast1,Correlation1,Energy1,Homogeneity1, Mean1, Standard_Deviation1, Entropy1, Variance1, Smoothness1, Kurtosis1, Skewness1, IDM1];


load svmStructColorRIRU_Final
load svmStructColorRUYE_Final
load svmStructColorRIYE_Final
 
result = zeros(1,3);
tmpResult = predict(svmStructColorRIRU,feat_color_disease);
if(tmpResult == 0)
    result(1,1) = result(1,1) + 1;
else
    result(1,2) = result(1,2) + 1;
end
tmpResult = predict(svmStructColorRUYE,feat_color_disease);
if(tmpResult == 1)
    result(1,2) = result(1,2) + 1;
else
    result(1,3) = result(1,3) + 1;
end
tmpResult = predict(svmStructColorRIYE,feat_color_disease);
if(tmpResult == 0)
    result(1,1) = result(1,1) + 1;
else
    result(1,3) = result(1,3) + 1;
end
 
[maxValue, maxIndex] = max(result);
switch maxIndex
    case 1
        finalResult = 'Ring Spot';
        pause(1)
        msgbox({'DISEASE AFFECTED:        Ring Spot';'';'';' CAUSING AGENT:        Leptosphaeria sacchari Phyllosticta .';'';'';' 	SYMPTOMS ARE';'';'';'1) Initial symptoms of ring spot are small, elongated, oval-shaped spots that are dark olivaceous green to reddish- brown with narrow yellow halos ';'';'';'2) Typical older symptoms of ring spot are large and elongated lesions (2.5 to 5 mm x 10 to 18 mm) with irregular outlines and red-brown margins ';'';'';'';});
            
    case 2
        finalResult = 'Sugar cane Rust';
        pause(1)
        msgbox({'DISEASE AFFECTED:       Sugar cane Rust';'';'';' CAUSING AGENT:        Puccinia melanocephala ';'';'';' 	SYMPTOMS ARE';'';'';'1) The earliest symptoms of common rust on the leaves are small, elongated yellowish spots which are visible on both the surfaces. These spots increase in size, mainly in length, and turn red-brown to brown in color. A narrow, pale yellow-green halo develops around the lesions.';'';'';' 2) Affected crops take on a distinctive rusty to brown color, as disease intensity builds up. The discoloration of these crops is readily seen from a distance.';'';'';'';});
        
    case 3
        finalResult = 'Yellow Spot';
        pause(1)
        msgbox({'DISEASE AFFECTED:        Yellow spot';'';'';' CAUSING AGENT:        Cercospora Koepkei';'';'';' 	SYMPTOMS ARE';'';'';'1) Presence of small, yellow coloured, irregularly shaped spots over the leaf surface. Density of spots is minimum in the lower surface, moderate in the middle and maximum towards the tip of the leaf.';'';'';'2) Spots coalesce at late stages and cause drying of leaves. Badly affected foliage looks reddish-brown when viewed from a distance.';'';'';'';});
        
end




% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I19=handles.ImgData6;
glcms = graycomatrix(I19);
% Derive Statistics from GLCM
stats = graycoprops(glcms,'Contrast Correlation Energy Homogeneity');
Contrast = stats.Contrast;
Correlation = stats.Correlation;
Energy = stats.Energy;
Homogeneity = stats.Homogeneity;
Mean = mean2(I19);
Standard_Deviation = std2(I19);
Entropy = entropy(I19);
RMS = mean2(rms(I19));
Variance = mean2(var(double(I19)));
a = sum(double(I19(:)));
Smoothness = 1-(1/(1+a));
Kurtosis = kurtosis(double(I19(:)));
Skewness = skewness(double(I19(:)));
% Inverse Difference Movement
m = size(I19,1);
n = size(I19,2);
in_diff = 0;
for i = 1:m
    for j = 1:n
        temp = I19(i,j)./(1+(i-j).^2);
        in_diff = in_diff+temp;
    end
end
IDM = double(in_diff);
feat_disease = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, RMS, Variance, Smoothness, Kurtosis, Skewness, IDM];

load svmStructRIRU
load svmStructRUYE
load svmStructRIYE
 
result = zeros(1,3);
tmpResult = predict(svmStructRIRU,feat_disease); 
if(tmpResult == 0)
    result(1,1) = result(1,1) + 1;
else
    result(1,2) = result(1,2) + 1;
end
tmpResult = predict(svmStructRUYE,feat_disease);
if(tmpResult == 1)
    result(1,2) = result(1,2) + 1;
else
    result(1,3) = result(1,3) + 1;
end
tmpResult = predict(svmStructRIYE,feat_disease);
if(tmpResult == 0)
    result(1,1) = result(1,1) + 1;
else
    result(1,3) = result(1,3) + 1;
end
 
[maxValue, maxIndex] = max(result);
switch maxIndex
    case 1
        finalResult = 'Ring Spot';
        pause(1)
        msgbox({'DISEASE AFFECTED:        Ring Spot';'';'';' CAUSING AGENT:        Leptosphaeria sacchari Phyllosticta .';'';'';' 	SYMPTOMS ARE';'';'';'1) Initial symptoms of ring spot are small, elongated, oval-shaped spots that are dark olivaceous green to reddish- brown with narrow yellow halos ';'';'';'2) Typical older symptoms of ring spot are large and elongated lesions (2.5 to 5 mm x 10 to 18 mm) with irregular outlines and red-brown margins ';'';'';'';});
            
    case 2
        finalResult = 'Sugar cane Rust';
        pause(1)
        msgbox({'DISEASE AFFECTED:       Sugar cane Rust';'';'';' CAUSING AGENT:        Puccinia melanocephala ';'';'';' 	SYMPTOMS ARE';'';'';'1) The earliest symptoms of common rust on the leaves are small, elongated yellowish spots which are visible on both the surfaces. These spots increase in size, mainly in length, and turn red-brown to brown in color. A narrow, pale yellow-green halo develops around the lesions.';'';'';' 2) Affected crops take on a distinctive rusty to brown color, as disease intensity builds up. The discoloration of these crops is readily seen from a distance.';'';'';'';});
        
    case 3
        finalResult = 'Yellow Spot';
        pause(1)
        msgbox({'DISEASE AFFECTED:        Yellow spot';'';'';' CAUSING AGENT:        Cercospora Koepkei';'';'';' 	SYMPTOMS ARE';'';'';'1) Presence of small, yellow coloured, irregularly shaped spots over the leaf surface. Density of spots is minimum in the lower surface, moderate in the middle and maximum towards the tip of the leaf.';'';'';'2) Spots coalesce at late stages and cause drying of leaves. Badly affected foliage looks reddish-brown when viewed from a distance.';'';'';'';});
        
end





function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
