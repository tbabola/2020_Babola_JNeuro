%%requires image processing toolbox, 
clear all; close all;
addpath(genpath('.\MATLAB Functions'));
[fn, dname] = uigetfile('*.czi'); %prompt for file to analyze
%open image

%%
imgdata = bfopen([dname fn]);
imgdata = imgdata{1,1}; %this is where the imagedata is stored
[fp,name,~] = fileparts([dname fn]); %store filename and path for saving later
oimg = []; index = 1;
for i = 1:1:size(imgdata,1) %because there are two channels in your image, only select GCaMP signal
    oimg(:,:,index) = imgdata{i,1};
    index = index + 1;
end
%%
%img = oimg(:,:,1:600); %alter this if drug condition;
%img = oimg(:,:,1200:1800); %MRS2500
img = oimg;
figure; imagesc(std(img,[],3));
if exist('IHCstruct')
    drug = 'MRS2500';
    drug = 'whole';
    points = IHCstruct.roiPoints;
    for i=1:size(points,1)
        drawpoint('Position',[points(i,1) points(i,2)]);
    end
else
    drug = [];
end
%%
if exist('IHCstruct')
    figure; imagesc(std(img,[],3));
    if exist('IHCstruct')
        for i=1:size(points,1)
            drawpoint('Position',[points(i,1) points(i,2)]);
        end
    end

    drug = 'MRS2500';
    oPoints = points;
end
%%
%redraw points


%% if fresh points
figure; imagesc(std(img,[],3));

points = []; roi = [];
while 1
    temp = drawpoint;
    if isempty(temp.Position)
        break
    else
        roi = [roi; temp];
        points = [points; temp.Position];
    end
end
oPoints = points;
%% reset points
% points = [];
% for i = 1:size(roi,1)
%     temp = roi(i);
%     points(i,:) = temp.Position;
% end
% 
% oPoints = points;
%%
movAvg = 5;
slope = 1;
points = oPoints;
points(:,2) = 512 - points(:,2) ;

for i=1:size(points,1)
    startPt = i - floor(movAvg/2);
    endPt = i + floor(movAvg/2);
    startPt;
    if startPt < 1
        startPt = 1;
    elseif endPt > size(points,1)
        endPt = size(points,1);
    end
    p = polyfit(points(startPt:endPt,1),points(startPt:endPt,2),1);
    slope = [slope; p(1)];
end

angles = atan(slope) * 180/pi;
angles = angles + 90;

%%
roisignal = zeros(size(img,3),size(points,1));
for i=1:size(points,1)
    h = drawellipse('Center',oPoints(i,:),'RotationAngle',angles(i),'SemiAxes',[3 3]);
    mask = h.createMask;
    I{i} = find(mask);
    
%     parfor j = 1:size(img,3)
%         temp = img(:,:,j);
%         roisignal(j,i) = mean(temp(I),'all');
%     end
end

for i=1:size(img,3)
    temp = img(:,:,i);
    for j=1:size(I,2)
        tempI = I{j};
        roisignal(i,j) = mean(temp(tempI),'all');
    end
end

%%

roisignalNorm = (roisignal - median(roisignal,1)) ./ median(roisignal,1);
roismooth = zeros(size(roisignalNorm));
for i=1:size(roisignalNorm,2)
    roismooth(:,i) = smooth(roisignalNorm(:,i));
end

medianroi = median(roismooth,1);
stdroi = std(roismooth,1);
thrroi = roismooth > medianroi + 2*stdroi;

%%
generateActivityMovieIHCs(Kalman_Stack_Filter(single(img)),thrroi,oPoints,angles,[dname '\' name '_activityMovie' drug],[500 25000])

%%
IHCstruct = struct();
IHCstruct.roisignal = roisignal;
IHCstruct.smoothroi = roismooth;
IHCstruct.thrroi = thrroi;
IHCstruct.angles = angles;
IHCstruct.roiPoints = oPoints;
save([dname '\' name '_IHCstruct' drug '.mat'],'IHCstruct');
