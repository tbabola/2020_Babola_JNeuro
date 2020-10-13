%% Snap25-T2A-GCaMP6s developmental video
%%% P1: 51-1250, P7:101-1301, P13: 3701-4900
file = '.\Data\Snap25_Dev_Stacked\Developmental_Snap25_Stacked.tif';
img = loadTif(file,32);

%%
leaderFrame = imread('.\Data\LeaderFrames\Developmental_Snap25_Stacked_leader.tif');
leaderFrame = imresize(leaderFrame,[NaN 480]);
imagesc(leaderFrame)
%%
name = ['.\Videos\SV08_Snap25_IC_development_HD'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 100;
end
v.FrameRate = 10;
open(v);
timernum = 0;
oimg = []; index = 1;
figure; imagesc(zeros(size(img,1,2)));  truesize; clims = [-10 300]; box off; axis off; axis off;
for i = -40:2:600
    if i <= 0
        imagesc(leaderFrame);
        axis off;
    else
        imagesc(img(:,:,i));
        colormap(gfb); axis off;
    end
    caxis(clims); 
    
    if i >= -40
        thandle = text(10, 15,'P1','FontSize',6,'Color','w','FontWeight','bold');
        thandle2 = text(10,162+15,'P7','FontSize',6,'Color','w','FontWeight','bold');
        thandle3 = text(10,162*2+20,'P13','FontSize',6,'Color','w','FontWeight','bold');
        if i >= 1
        thandle4 = text(5,480-5,'\it Snap25-T2A-GCaMP6s','FontSize',5,'Color','w','FontWeight','bold','HorizontalAlignment','left');
        end
    end
    
    if i > 1 & mod(i,10) == 0
        timernum = i/10;
    end
    if i > 1
        timer = text(480-10,485-15,[num2str(timernum) ' s'],'FontSize',6,'color','w','FontWeight','bold','HorizontalAlignment','right');
    end
    
    if mod(i,100) == 0
        disp(i);
    end
    writeVideo(v,getframe); 
end

close(v);

%% Snap25-T2A-GCaMP6s a9 KO/GOF/control
clear all; close all;
file = '.\Data\a9_Snap25_Stacked\a9_Snap25_Stacked_2.tif';
img = loadTif(file,32);

%%
leaderFrame = imread('.\Data\LeaderFrames\Snap25_a9_leader_2.tif');
imagesc(leaderFrame)

%%
name = ['.\Videos\SV09_a9KO_GOF_HD'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 100;
end
v.FrameRate = 10;
open(v);
timernum = 0;
oimg = []; index = 1;
figure; imagesc(zeros(size(img,1,2)));  truesize; clims = [-10 450]; box off; axis off; axis off;
for i = -40:2:600
    if i <= 0
        imagesc(leaderFrame);
        axis off;
    else
        imagesc(img(:,:,i));
        colormap(gfb); axis off;
    end
    caxis(clims); 
    
    if i >= -40
        thandle = text(5, 145-10,'\alpha\it9 \rm\bfKO','FontSize',8,'Color','w','FontWeight','bold');
        thandle2 = text(5,145*2,'Control','FontSize',8,'Color','w','FontWeight','bold');
        thandle3 = text(5,460-15,'\alpha\it9 \rm\bfGOF','FontSize',8,'Color','w','FontWeight','bold');
    end
    
    if i > 1 & mod(i,10) == 0
        timernum = i/10;
    end
    if i > 1
        timer = text(480-10,465-15,[num2str(timernum) ' s'],'FontSize',8,'color','w','FontWeight','bold','HorizontalAlignment','right');
    end
    
    if mod(i,100) == 0
        disp(i);
    end
    writeVideo(v,getframe); 
end

close(v);

