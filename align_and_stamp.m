clear all
%Program to align 2P images using their coordinates

%Images must first be made into max projections using imageJ to output to a
%sequence of images using the numbering '000,001',etc. This is because
%matlab cannot read the raw TIFFs from the 2P system

%Coordinates must be copied to a separate TXT file in order of the image
%sequence, i.e. the image no. will not be read
%ALSO, coordinates must be in the format (NUMBER,NUMBER) with no other
%brackets in the text file.

zoom_factor=150; 
%ATTENTION: This is the number of microns in one dimension of a single
%image. For z2, this is 150um, for z3 this is 100um.

%----------------------------------------------------------------------------%
%STEP ONE: Get filename for importing a series of images & directory
%e.g. Pup1_LH2_Session_a_001_AVG.jpg to =imread(''Pup1_LH3_b-',num2str(i,'%03i'),'_AVG.jpg'');
%must be in the format of 000,001,002, etc.)
[name,pathstr]=uigetfile('*tif','Select first image in the series, e.g. Pup1_LH2_a_001.jpg'); %user selects first file

cd(pathstr); %set path as the directory with the first image in it (optional)

filename1_idx=strfind(name,'00');%finds the number in the filename
filename1=name(1:filename1_idx-1);%finds the first part of the filename, up to the number
filename2=name(filename1_idx+3:end);%finds the second part of the filename, after the number

%----------------------------------------------------------------------------%
%STEP TWO: Grab coordinates from within brackets of a text file
document=uigetfile('*.txt', 'Select the file with the coordinates in brackets, in order of image seq.'); %gets the file
fid=fopen(document);
f=fscanf(fid,'%s'); %turns all text into a single string
index_bracket= strfind(f,'('); %finds ( in the string
index_bracket2= strfind(f,')'); %finds ) in the string
[n,m]=size(index_bracket); % m is now how many values there are in the whole document
[n2,m2]=size(f); %m2 is now how many characters there are in the string, f

index_bracket=index_bracket'; %transpose to a list
index_bracket2=index_bracket2'; %transpose to a list

index_bracket=index_bracket+1; %move place to first number inside brackets
index_bracket2=index_bracket2-1;%move place to last number inside brackets

fx=[0 0]; %set variables
fy=[0 0]; %set variables

for i=1:m 
value1=str2num(f(index_bracket(i):index_bracket2(i))); %writes coordinates into var:Value1
co_x=value1(:,1);
co_y=value1(:,2);
fx=[fx,co_x]; %writes separate coordinates to fx or fy
fy=[fy,co_y];
end

fx=fx(3:m+2);%removes pre-allocation
fy=fy(3:m+2);%removes pre-allocation

coords_x=fx';%transpose to list and set as var: coords_x/y
coords_y=fy';

%----------------------------------------------------------------------------%

%STEP THREE: Convert coordinates to pixels
%i.e. coordinate width & height divided by 512
coord_convert=512/zoom_factor; %this is the number of pixels in one coordinate value, i.e. 512/one length of image in coords. N.b. 150 is only for zoom2x
coords_x=coords_x*coord_convert;
coords_y=coords_y*coord_convert;%converts values

%   b. Take minimum and maximum XY values to make complete matrix

x_min=min(coords_x);
y_min=min(coords_y);
%   c. Zero data to 0,0 and above

coords_x=coords_x-x_min;
coords_y=coords_y-y_min;
x_max=max(coords_x);
y_max=max(coords_y);
coords_x=ceil(coords_x); %rounds all numbers up to full integers
coords_y=ceil(coords_y);

%y coordinates need to be 'reversed' for axons going from top right to
%bottom left. This is MAX - the value

coords_y=max(coords_y)-coords_y; %IMPORTANT: This reverses the y-values, for some reason this seems to be needed to order the coordinates properly. Seems to be the same process for any axon...
%coords_x=max(coords_x)-coords_x;

%----------------------------------------------------------------------------%

%STEP FOUR: Make the total map, pre-allocated with zeros

total_map=zeros(ceil(y_max),ceil(x_max)); %makes zeroes matrix of X and Y maximum values, ceiling is because you can't have '0.5' of a pixel

total_map=uint8(total_map); %convert to uint8 format, this should be whatever form the original images are in

%----------------------------------------------------------------------------%

%STEP FIVE: Import images
x_number=size(coords_x);
x_number=x_number(:,1); %reads number of images

for i=1:x_number %change number to starting number of images, as selected in line 21, default =1
   str1=['''',filename1,num2str(i,'%03i'),filename2,'''']; %creates the full filename with 'i' being the number
   eval(['image=imread(',str1,')',';',''])
   
   
        %-----------INSERTION OF NAME-STAMPING CODE------------%
   %Makes white image the same size as the oringinal
hf = figure('color','white','units','normalized','position',[.1 .1 .8 .8]);
image(ones(size(image))); 
set(gca,'units','pixels','position',[5 5 size(image,2)-1 size(image,1)-1],'visible','off')


% Text at arbitrary position 
text('units','pixels','position',[12 200],'rotation',90,'fontsize',15,'string',num2str(i,'%03i')) 

% Capture the text image 
% Note that the size will have changed by about 1 pixel 
tim = getframe(gca); 
close(hf) 

% Extract the cdata
tim2 = tim.cdata;

% Make a mask with the negative of the text 
tmask = tim2==0; 

% Place white text 
% Replace mask pixels with UINT8 max 
tmask=tmask(:,:,1);
tim2=tmask(:,:,1);
image(tmask) = uint8(255); 

            %---------------------------------------------------%
           

   eval(['img_',num2str(i,'%03i'),'=image;']); %loads images named '00#.jpg', where # is 1:the number of images. Saves these as 'img_00#' sequentially
end

%----------------------------------------------------------------------------%

%STEP SIX: Copies images to the coordinates on the total map

j=1;
for i=1:x_number
  eval(['place_img=img_',num2str(i,'%03i'),';'])
   total_map(((coords_y(j))+1):((coords_y(j))+512),((coords_x(j))+1):((coords_x(j))+512))=place_img;
   j=j+1;
end

%----------------------------------------------------------------------------%
%Show final image and write to file
imshow(total_map)

if zoom_factor==150
    savezoom='complete_z2.tif';
elseif zoom_factor==100
    savezoom='complete_z2.tif';
else 
     savezoom=['complete_', num2str(zoom_factor), 'z2.tif'];
     end
imwrite(total_map,savezoom,'TIF')




