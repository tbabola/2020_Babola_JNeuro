%%
clear all; close all;
fileList= loadFileList('.\Data\Snap25_MRS2500\*_*.tif');
baseline = loadTif(fileList{1},16);
MRS= loadTif(fileList{2},16);
%%
writeTif(int16(mean([baseline(:,:,10:19) MRS(:,:,1:10)],3)),'.\Data\Snap25_MRS2500\firstFrame.tif',16)
%%
leaderFrame = imread('.\Data\LeaderFrames\Snap25_MRS2500_leader.tif');
figure; imagesc(leaderFrame);
%%
name = ['.\Videos\SV07_SnapMRS2500_comp'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 55;
end
v.FrameRate = 10;
open(v);

green = zeros(255,3);
green(:,2) = [1:255]';
green = green/255;

drugApp = 300; drugOff = drugApp + 60*15;
figure; imagesc(zeros(240,480)); colormap(gray); truesize; box off; hold on; axis off;
for i = -10:1:300
    if i <= 0
        imagesc(imresize(leaderFrame,[NaN 480])); 
        colormap(green);
    else
        imagesc(imresize([baseline(:,:,i) MRS(:,:,i)],[NaN 480])); caxis([-200 3000])
    end

    if i == -10
        thandle0 = text(5,225,'P0 Base','FontSize',12,'color','w','FontWeight','bold')
        thandle = text(240-5, 10,'Baseline','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','right');
        thandle2 = text(480-5, 10,'+MRS2500','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','right');
    else
        set(thandle0,'Layer', 'front');
        set(thandle ,'Layer', 'front');
        set(thandle2 ,'Layer', 'front');
        if exist('timer')
            delete(timer)
        end
    end
    if i > 1
        timer = text(240,225,[num2str(i) ' s'],'FontSize',12,'color','w','FontWeight','bold','HorizontalAlignment','center');
    end
    
    if mod(i,100) == 0
        disp(i);
    end
    writeVideo(v,getframe); 
end

close(v);
