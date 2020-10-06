%%MRS2500 Video
%video to use
clear all; close all;
fileList = loadFileList('.\Data\Pax2Cre\*45*.czi');
imgdata = bfopen(fileList{1});
imgdata = imgdata{1,1}; %this is where the imagedata is stored
%%
img = []; index = 1;
parfor i = 1:2400
    img(:,:,i) = imgdata{i,1};
end
clear imgdata;
img = imresize(img,[NaN 480]);
img = Kalman_Stack_Filter(img);
%%
%writeTif(single(mean(img(:,:,1:5),3)),'.\VidProcessing\firstFrameMRSandAnalysis.tif',32);
%%
leaderFrame = imread('.\Data\LeaderFrames\Pax2_MRS_leader.tif');
%leaderFrame = leaderFrame(:,:,1:3);
%leaderFrame
%%
name = ['.\Videos\SV03_Pax2_MRS2500_comp'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 34;
end
v.FrameRate = 10;
open(v);

green = zeros(255,3);
green(:,2) = [1:255]';
green = green/255;


oimg = []; index = 1; drugApp = 825;
timernum = '0';
figure; imagesc(zeros(480,480)); truesize; clims = [500 30000]; box off; hold on;
for i = 560:2:1200
    if i < 600
        imagesc(leaderFrame);
    elseif i >= 600
        hold off;
        imagesc(img(:,:,i)); colormap(green); caxis(clims); 
        hold on;
        if i >= drugApp
            thandle = text(10, 40,'+MRS2500','FontSize',16,'Color','w','FontWeight','bold');
        end
        
        if mod(i,2) == 0 
          timernum = num2str((i-600)/2);
        end
        timer = text(480-10,480-15,[timernum ' s'],'FontSize',16,'color','w','FontWeight','bold','HorizontalAlignment','right');

        if mod(i,100) == 0
            disp(i);
        end
    end
     
    if i >= 560
        thandle0 = text(10,15,'P0\it Pax2-Cre;GCaMP3','FontSize',16,'color','w','FontWeight','bold');
    end
    writeVideo(v,getframe); 
end

close(v);


%% might as well use same video for analysis video
name = ['.\Videos\SV02_ISCIHCanalysis_comp'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 29;
end
v.FrameRate = 10;
open(v);

green = zeros(255,3);
green(:,2) = [1:255]';
green = green/255;

fileList = loadFileList('.\Data\Pax2Cre\*45*.mat')
%%%img,labels,positiveIndices,IHClabelroi,IHClocation,name,clims,h
h = load(fileList{1});
IHCstruct = h.IHCstruct;
h = load(fileList{2});
ISCstruct = h.ISCstruct;
labels = ISCstruct.labelRoi;
positiveIndices = ISCstruct.posIndices*480/512;
IHClabelroi = IHCstruct.thrroi;
IHClocation = IHCstruct.IHCcenters*480/512;
timernum = 0;
figure; imagesc(zeros(480,480)); truesize; clims = [500 30000]; box off; hold on;
for i = -20:1:300
     if i <= 0
        imagesc(leaderFrame);
        if exist('timer')
            delete(timer)
        end
    else
        indx_green = find(labels(i,:)>=1);
            imagesc(img(:,:,i)); colormap(green); caxis(clims); hold on;
            if ~isempty(indx_green)
                for j = 1:size(indx_green,2)
                     plot(positiveIndices(indx_green(j),1:2:end),positiveIndices(indx_green(j),2:2:end),'Color','w');
                end
            end

            %%IHCs
            indx_green = find(IHClabelroi(i,:)>=1);
            if ~isempty(indx_green)
                for j = 1:size(indx_green,2)
                    temp = indx_green(j);
                    drawellipse('Position',IHClocation(temp,:),'SemiAxes',[1 1],'Color','w','InteractionsAllowed','none');
                end
            end          
            

            if mod(i,2) == 0 
                timernum = num2str(i/2);
            end

            timer = text(480-10,480-15,[timernum ' s'],'FontSize',16,'color','w','FontWeight','bold','HorizontalAlignment','right');

            drawnow;
            
            hold off;
            if mod(i,100) == 0
                disp(i);
            end
     end

    thandle0 = text(10,15,'P0\it Pax2-Cre;GCaMP3','FontSize',16,'color','w','FontWeight','bold');
    writeVideo(v,getframe); 
end   

close(v)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%apex/base video
%Image 41 and Image 39, same as figure
clear all; close all;
fileList = loadFileList('.\Data\Pax2Cre\*35*.czi');
imgdata = bfopen(fileList{1});
imgdata = imgdata{1,1}; %this is where the imagedata is stored
img = []; index = 1;
parfor i = 1:600
    img(:,:,i) = imgdata{i,1};
end
clear imgdata;
img = Kalman_Stack_Filter(img);
%%
fileList = loadFileList('.\Data\Pax2Cre\*41*.czi');
imgdata = bfopen(fileList{1});
imgdata = imgdata{1,1}; %this is where the imagedata is stored
imgapex = []; index = 1;
parfor i = 1:600
    imgapex(:,:,i) = imgdata{i,1};
end
clear imgdata;
imgapex = Kalman_Stack_Filter(imgapex);
%%
bigImg = [imgapex img];
clear img imgapex
bigImg = imresize(bigImg,[NaN 480]);
%%
%writeTif(single(mean(bigImg(:,:,1:5),3)),'.\VidProcessing\firstFrameBaseApex.tif',32)
%%
leaderFrame = imread('.\Data\LeaderFrames\Pax2_Base_Apex_leader.tif');
%%
name = ['.\Videos\SV04_ApexvsBase_comp'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 37;
end
v.FrameRate = 10;
open(v);

green = zeros(255,3);
green(:,2) = [1:255]';
green = green/255;

oimg = []; index = 1;
timernum = 0;
figure; imagesc(zeros(size(bigImg,1,2))); colormap(green); truesize; clims = [500 30000]; box off; hold on; axis off;
for i = -20:1:600
    if i < 1
        imagesc(leaderFrame)
    else
        hold off; imagesc(bigImg(:,:,i)); caxis(clims); box off; axis off; hold on;
        if mod(i,2) == 0 
           timernum = num2str(i/2);
        end

        timer = text(480-5,240-15,[timernum ' s'],'FontSize',12,'color','w','FontWeight','bold','HorizontalAlignment','right');

        if mod(i,100) == 0
            disp(i);
        end
    end
    thandle = text(10, 15,'P0 Apex','FontSize',12,'Color','w','FontWeight','bold');
    thandle2 = text(240+10,15,'P0 Base','FontSize',12,'Color','w','FontWeight','bold');
    writeVideo(v,getframe); 
end

close(v);


