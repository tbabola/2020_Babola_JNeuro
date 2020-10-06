%%Crenation video
fileList= loadFileList('.\Data\Crenations\*.tif');

E16img = loadTif(fileList{4},8);
P0img = loadTif(fileList{1},8);
P7img = loadTif(fileList{2},8);
P11img = loadTif(fileList{3},8);

bigImg = [E16img P0img; P7img P11img+20];
clear E16img P0img P7img P11img
bigImg = imresize(bigImg, [NaN 480]);
%%
writeTif([E16img(:,:,1) P0img(:,:,1); P7img(:,:,1) P11img(:,:,1)+20],'.\Data\ProcessedVidsForFigure\firstFrame.tif',8)
%%
leaderFrame = loadTif('.\Data\LeaderFrames\Crenation_leader.tif');
leaderFrame = imresize(leaderFrame,[NaN 480]);
imagesc(leaderFrame)
%%
name = ['.\Videos\SV01_Crenations_comp'];
ext = '.mp4';
fileOutput = [name ext];
if strcmp(ext,'.avi')
    v = VideoWriter(fileOutput,'Uncompressed AVI');
elseif strcmp(ext,'.mp4')
    v = VideoWriter(fileOutput,'MPEG-4');
    v.Quality = 8;
end
v.FrameRate = 20;
open(v);

oimg = []; index = 1;
drugApp = 300; drugOff = drugApp + 60*15;
figure; imagesc(zeros(size(bigImg,1,2))); colormap(gray); truesize; clims = [0 200]; box off; hold on; axis off;
for i = -40:2:2100
    if i <= 0
        imagesc(leaderFrame);
    else
        imagesc(bigImg(:,:,i))
    end
    caxis(clims); 
    
    if i == -40
        thandle = text(5, 10,'E16','FontSize',12,'Color','w','FontWeight','bold');
        thandle2 = text(240+5,10,'P0','FontSize',12,'Color','w','FontWeight','bold');
        thandle3 = text(5,139+10,'P7','FontSize',12,'Color','w','FontWeight','bold');
        thandle4 = text(240+5,139+10,'P11','FontSize',12,'Color','w','FontWeight','bold');
    else
        set(thandle ,'Layer', 'front')
        set(thandle2 ,'Layer', 'front')
        set(thandle3 ,'Layer', 'front')
        set(thandle4,'Layer', 'front')
    end
    
    if i == drugApp
        thandle5 = text(240 - 5, 10,'+MRS2500','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
        thandle6 = text(480 - 5,10,'+MRS2500','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
        thandle7 = text(240 - 5,139+10,'+MRS2500','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
        thandle8 = text(480 - 5,139+10,'+MRS2500','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
    elseif i > drugApp & i < drugOff
        set(thandle5 ,'Layer', 'front')
        set(thandle6 ,'Layer', 'front')
        set(thandle7 ,'Layer', 'front')
        set(thandle8,'Layer', 'front')
    elseif i == drugOff
        delete(thandle5);
        delete(thandle6);
        delete(thandle7);
        delete(thandle8);
        thandle5 = text(240 - 5, 10,'Washout','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
        thandle6 = text(480 - 5,10,'Washout','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
        thandle7 = text(240 - 5,139+10,'Washout','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
        thandle8 = text(480 - 5,139+10,'Washout','FontSize',12,'Color','w','FontWeight','bold','HorizontalAlignment','Right');
    elseif i > drugOff
        set(thandle5 ,'Layer', 'front')
        set(thandle6 ,'Layer', 'front')
        set(thandle7 ,'Layer', 'front')
        set(thandle8,'Layer', 'front')
        
    end
    
    if i > 1
        timer = text(240,270,[num2str(i) ' s'],'FontSize',12,'color','w','FontWeight','bold','HorizontalAlignment','center');
    end
    
    if mod(i,100) == 0
        disp(i);
    end
    writeVideo(v,getframe); 
end

close(v);
