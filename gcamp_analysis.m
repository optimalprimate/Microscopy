%Program to analyse TIF files from t-series from GCAMP6 imaging
%Aims: 
    %import files, to time stack, mask pixels based on activity/excl.
    %over-sat
    %plot intensity over time
    %FFT pattern recognition
    %quantitative - activity value for neuron
    %potential to measure red channel
    %export to video
    
    
    
clear all
close all

%set frames per second
%fps=input('what FPS was used?   ');
fps=3.4;
fps=fps/2;
disp('Select Images...')

%stage1 - import of TIFs
[name,pathstr]=uigetfile('*tif','Select the first TIF in the folder'); %user selects first file

cd(pathstr); %set path as the directory with the first img in it (optional)
filelist=dir('*.tif'); %lists TIF files in dir

filename1_idx=strfind(name,'Cycle');%finds the text Cycle in filename
filename1=name(1:filename1_idx+4); %loads first part up to the cycle number
filename2=name(filename1_idx+10:end); %loads second part of file name after cycle number

%----------------------------------------------------------------------------%
%----------------------------------------------------------------------------%

h = waitbar(0,'Loading Images...');
%load imgs - creates a stack
for i=1:length(filelist)
    str1=['''',filename1,num2str(i,'%05i'),filename2,'''']; %creates the full filename with 'i' being the number
    eval(['img(:,:,',num2str(i),')=imread(',str1,')',';',''])
    
    %load notes
perc = 100/length(filelist)*i;
waitbar(perc/100,h,sprintf('%d%% complete...',ceil(perc)))
    
end

close(h)


%threshold pixels to remove background (based on average intensity px)
stack_size=size(img);
disp('Averaging...')
%takes the average of each pixel
for j=1:stack_size(2) %y-axis
    for i=1:stack_size(1) %x-axis of each img
        img2=mean2(img(i,j,:));
        mean_img(i,j)=img2(1,:);
        end
    
end
disp('Thresholding...')
%replaces all pixels with a mean over with 0
img_thresholded=img;
for j=1:stack_size(2)
    for i=1:stack_size(1)
        if mean_img(i,j)<25
            img_thresholded(i,j,:)=0;
        else
            
        end
    end
end
disp('Select Rough Crop')
%display average image for crop and video to check activity
img_forplay = uint8(255*mat2gray(img));
implay(img_forplay,50)

figure('name','Select box around axon of interest (background will be removed)'); 
[mask,mask_coords]=imcrop(mean_img,[]); %mask_coords is now the coordinates for the corners of the ROI
 close all
 
 disp('Masking Images...')
%crop stack to mask + average of each img
for i=1:length(filelist)
 img_cropped(:,:,i)=imcrop(img_thresholded(:,:,i),mask_coords);
 z=img_cropped(:,:,i); %take single image for...
 z(z==0)=[]; %remove threhsolded 0s from single image for mean
 means(i)=mean2(z);
 z=double(z);
 variances(i)=var(z);
 %means(i)=mean2(img_thresholded(:,:,i));
end


disp('Displaying Analysis.')

x=1000/fps; %each frame in s


%plot


figresults=figure('name',name(1:end-30));

xaxis=1:x:(x*length(filelist));
xaxis=xaxis/1000;
subplot(2,2,1)
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
plot(xaxis,means)
title('Mean intensity over time')
xlabel('Time (s)')
subplot (2,2,2)
plot(xaxis,smooth(means, 8, 'moving'), 'r' )
title('Moving average of mean intensity over time')
xlabel('Time (s)')

%textbox to display title
mTextBox = uicontrol('style','text');
figPosition = get(figresults,'Position');
y=figPosition(4)-30;
set(mTextBox,'String',name(1:end-30),'Position',[200 y 400 20])
%mTextBoxPosition = get(mTextBox,'Position');
%set(mTextBox,'Units','characters')


%img_cropped_8=uint8(img_cropped);
img_cropped8 = uint8(255*mat2gray(img_cropped));
implay(img_cropped8)

%show variances
%subplot(2,2,3)
%plot(xaxis,variances)
%title('Variance over time')
%xlabel('Time (s)')

%show mini montage
j=1;
cropsize=size(img_cropped);
imgresh=reshape(img_cropped,cropsize(1),cropsize(2), 1, cropsize(3));
for i = 1:20:length(filelist)
    imgmont(:,:,:,j)=imgresh(:,:,:,i);
    j=j+1;
end
imgmont = uint8(255*mat2gray(imgmont));
subplot(2,2,3)
montage(imgmont)
title('Montage of every 20th image')





%FFT analysis

x=x/1000;
Y = fft(means);
Y(1) = [];
n = length(Y);
power = abs(Y(1:floor(n/2))).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;
period=x./freq;
subplot(2,2,4)
plot(period(1:20),power(1:20))
xlabel('Frequency')
title('Fourier Transform')
hold on;
index = find(power == max(power));
if period(index) > (length(filelist)/fps)-5
    sorted=sort(power);
    find_index=sorted([end-1]);
    index = find(power==find_index);
end

maxfreq=period(index);
mainPeriodStr = num2str(period(index));
plot(period(index),power(index),'r.', 'MarkerSize',25);
text(period(index)+2,power(index),['Most frequent event every ',mainPeriodStr,' seconds']);
hold off;

%write analysis
cd ..
if exist('Analysis') ~=7
mkdir ('Analysis')
else
end

cd Analysis

filenamexls=[filename1,'.xls'];
smoothmeans=smooth(means, 8, 'moving');
xlsdataheaders={'Time', 'Average', 'Smooth Average', 'FFT Peak','FPS'};
for i=1:length(means)
      xlsdata_numbers(:,i)=[xaxis(i),means(i),smoothmeans(i),maxfreq,fps];

end

xlsdata_numbers= num2cell(xlsdata_numbers');
xls_export=[xlsdataheaders;xlsdata_numbers];

xlswrite(filenamexls,xls_export)

    








