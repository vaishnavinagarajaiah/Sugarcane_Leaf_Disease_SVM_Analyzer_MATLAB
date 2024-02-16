clc;
%clear all;
%ii=1;
%jj=1;
%while jj<12

[filename, pathname] = uigetfile({'*.*';'*.bmp';'*.jpg';'*.gif'}, 'Pick a Leaf Image File');
%image=int2str(jj);
%image=strcat(int2str(jj),'.jpg');
%I=imread(image);
I = imread([pathname,filename]);

%figure,imshow(I),title('original');
I=imresize(I,[255 255]);

%lab_he=rgb2hsv(I);
cform = makecform('srgb2lab');
lab_he = applycform(I,cform);
figure,imshow(lab_he),title('L*a*b');
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);


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

figure,imshow(segmented_images{1}), title('objects in cluster 1');
figure,imshow(segmented_images{2}), title('objects in cluster 2');
figure,imshow(segmented_images{3}), title('objects in cluster 3');

%mean_cluster_value = mean(cluster_center,2);
%[tmp, idx] = sort(mean_cluster_value);
%blue_cluster_num = idx(1);

%L = lab_he(:,:,1);
%%blue_idx = find(pixel_labels == blue_cluster_num);
%L_blue = L(blue_idx);
%is_light_blue = im2bw(L_blue,graythresh(L_blue));

%nuclei_labels = repmat(uint8(0),[nrows ncols]);
%nuclei_labels(blue_idx(is_light_blue==false)) = 1;
%nuclei_labels = repmat(nuclei_labels,[1 1 3]);
%blue_nuclei = I;
%blue_nuclei(nuclei_labels ~= 1) = 0;
%figure,imshow(blue_nuclei), title('blue nuclei');


%Res1=imsubtract(segmented_images{1},segmented_images{2});
%Res2=imsubtract(segmented_images{2},segmented_images{3});
%Res3=imsubtract(segmented_images{1},segmented_images{3});
%Res4=imsubtract(segmented_images{2},segmented_images{1});
%Res5=imsubtract(segmented_images{3},segmented_images{2});
%Res6=imsubtract(segmented_images{3},segmented_images{1});
%figure,subplot(2,3,1), imshow(Res1);
%subplot(2,3,2), imshow(Res2);
%subplot(2,3,3),imshow(Res3);
%subplot(2,3,4), imshow(Res4);
%subplot(2,3,5), imshow(Res5);
%subplot(2,3,6),imshow(Res6);
%add=imadd(Res1,Res2);
%figure, imshow(add);
%im=add;
prompt='enter value from 1 to 3';
x=input(prompt);
%k=1;
%x=1;
%while k<4

g1=rgb2gray(segmented_images{k});
e1=edge(g1,'sobel');
figure,imshow(e1);
I19=e1;
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
feat_color_disease = [Contrast,Correlation,Energy,Homogeneity, Mean, Standard_Deviation, Entropy, Variance, Smoothness, Kurtosis, Skewness, IDM];
%feat_color_variable(ii,1)=Contrast;
%feat_color_variable(ii,2)=Correlation;
%feat_color_variable(ii,3)=Energy;
%feat_color_variable(ii,4)=Homogeneity;
%feat_color_variable(ii,5)=Mean;
%feat_color_variable(ii,6)=Standard_Deviation;
%feat_color_variable(ii,7)=Entropy;
%feat_variable(ii,8)=RMS;
%feat_color_variable(ii,8)=Variance;
%feat_color_variable(ii,9)=Smoothness;
%feat_color_variable(ii,10)=Kurtosis;
%feat_color_variable(ii,11)=Skewness;
%feat_color_variable(ii,12)=IDM;
%k=k+1;
%ii=ii+1;
%end
%jj=jj+1;
%end

%for i=1:120
 %   if i<40
  %      FeatureColorLabel(i,:)=0;
   % end 
    %if i>39 && i<88
            
     %   FeatureColorLabel(i,1)=1;
    %end
    
    %if i>87
     %   FeatureColorLabel(i,1)=2;
    %end
        
%end
%save('FeatureColorLableTrain.mat','FeatureColorLabel');
load svmStructColorRIRU
load svmStructColorRUYE
load svmStructColorRIYE
 
result = zeros(1,3);
tmpResult = svmclassify(svmStructColorRIRU,feat_color_disease);
if(tmpResult == 0)
    result(1,1) = result(1,1) + 1;
else
    result(1,2) = result(1,2) + 1;
end
tmpResult = svmclassify(svmStructColorRUYE,feat_color_disease);
if(tmpResult == 1)
    result(1,2) = result(1,2) + 1;
else
    result(1,3) = result(1,3) + 1;
end
tmpResult = svmclassify(svmStructColorRIYE,feat_color_disease);
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
