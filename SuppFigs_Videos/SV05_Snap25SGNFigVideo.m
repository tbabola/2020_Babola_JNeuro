%% Snap25 SGN Supplemental Video
fileList= loadFileList('.\Data\Snap25_SGN_Vids\*_*.tif');
E16apex = loadTif(fileList{1},16);
E16base = loadTif(fileList{2},16);
P0apex = loadTif(fileList{3},16);
P0base = loadTif(fileList{4},16);
%%
writeTif(int16(mean([E16apex(:,:,1:10) P0apex(:,:,1:10); E16base(:,:,1:10) P0base(:,:,1:10)],3)),'.\Data\Snap25_SGN_Vids\firstFrame.tif',16)
%%
leaderFrame = imread('.\Data\LeaderFrames\Snap25_SGN_Stacked_leader.tif');
imagesc(leaderFrame)
leaderFrame(:,799:800,:) = 0;
leaderFrame(799:800,:,:) = 0;
%%
name = ['.\Videos\SV05_SnapDev_comp'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 35;
end
v.FrameRate = 10;
open(v);

green = zeros(255,3);
green(:,2) = [1:255]';
green = green/255;

drugApp = 300; drugOff = drugApp + 60*15;
figure; imagesc(zeros(size(imresize(leaderFrame,[NaN 480])))); colormap(gray); truesize; box off; hold on; axis off;
for i = -10:1:300%:size(E16img,3)
    if i <= 0
        imagesc(imresize(leaderFrame,[NaN 480])); 
        colormap(green);
    else
        imagesc(imresize([E16apex(:,:,i) P0apex(:,:,i); E16base(:,:,i) P0base(:,:,i)],[NaN 480])); caxis([-500 11000])
    end

    if i == -10
        thandle0 = text(5,10,'\it Snap25-GCaMP6s','FontSize',12,'color','w','FontWeight','bold')
        thandle = text(240-5, 10,'E16 apex','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','right');
        thandle2 = text(480-5,10,'P0 apex','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','right');
        thandle3 = text(240-5,240+10,'E16 base','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','right');
        thandle4 = text(480-5,240+10,'P0 base','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','right');
    else
        set(thandle0,'Layer', 'front');
        set(thandle ,'Layer', 'front');
        set(thandle2 ,'Layer', 'front');
        set(thandle3 ,'Layer', 'front');
        set(thandle4,'Layer', 'front');
        if exist('timer')
            delete(timer)
        end
    end
    if i > 1
        timer = text(240,465,[num2str(i) ' s'],'FontSize',12,'color','w','FontWeight','bold','HorizontalAlignment','center');
    end
    
    if mod(i,100) == 0
        disp(i);
    end
    writeVideo(v,getframe); 
end

close(v);