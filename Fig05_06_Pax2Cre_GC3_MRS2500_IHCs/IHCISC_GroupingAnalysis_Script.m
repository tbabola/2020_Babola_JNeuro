%load image
clear all; close all;
addpath(genpath('..\MATLAB Functions'));
[fn, dname] = uigetfile('*.czi'); %prompt for file to analyze
%open image
%%
imgdata = bfopen([dname fn]);
imgdata = imgdata{1,1}; %this is where the imagedata is stored
[fp,name,~] = fileparts([dname fn]); %store filename and path for saving later
img = []; index = 1;
for i = 1:600%1:size(imgdata,1) 
    img(:,:,index) = imgdata{i,1};
    index = index + 1;
end
[fp,name,~] = fileparts([dname fn]);
%%
figure(1); imagesc(mean(img,3));
img = bleachCorrect(img,2);
x = input('Rotate? (degrees)');
    ISCstruct.rotateDegrees = x;
    if x ~=0
        img = imrotate(img,x);
        imagesc(mean(img,3)); truesize;
    end
%%
ISCstruct = getISCevents(img,fp,name,0);

%%
close all;
IHCstruct = getIHCeventsAuto(img,fp,name,0);

%%
%now the real work begins
ISClabels = ISCstruct.labelRoi;
ISCrois = ISCstruct.rois;
IHCthr = IHCstruct.thrroi;
IHCrois = IHCstruct.roisignalNorm;
numEvents = max(ISClabels,[],'all');
timesOfEvents = [[ISCstruct.event(:).timeStart]' [ISCstruct.event(:).timeEnd]'];

%sort rois on x location
[~,newOrder] = sort(IHCstruct.bottomPos(:,1));
IHCthr = IHCthr(:,newOrder);
IHClabels = int16(IHCthr);
IHClabels(IHClabels == 1) = 999;
IHClocs = IHCstruct.bottomPos(newOrder,:);
figure; imagesc([IHCthr*10 bwlabel(imgaussfilt(single(IHCthr'),[2,0.5])>0.25)']);
rawIHCgrouplabels = bwlabel(IHCthr',4);
rawIHCgrouplabels = rawIHCgrouplabels';
IHCgrouplabels = bwlabel(imgaussfilt(single(IHCthr'),[2,0.5])>0.25);
IHCgrouplabels = IHCgrouplabels';
numLabels = max(IHCgrouplabels,[],'all');
h = figure(3);
%process the large IHC events
numLabelsOrig = numLabels;
for i = 1:numLabelsOrig
    [r,c] = find(IHCgrouplabels == i);
    activeIHCs = unique(c);
    if size(activeIHCs,1) > 0.6 * size(IHClabels,2) | (max(r) - min(r)) > 10
        time = [min(r) max(r)];
        tempIHClabels = IHCgrouplabels == i;
        generateISCIHCActivityMovie(img(:,:,time(1):time(2)),ISClabels(time(1):time(2),:),ISCstruct.posIndices,tempIHClabels(time(1):time(2),:),IHCstruct.IHCcenters(newOrder,:),[],[0 max(img,[],'all')],h)
        satisfied = 0;
        while ~satisfied
            x = input('How many events do you see? replay (0) ');
            if x == 0
                generateISCIHCActivityMovie(img(:,:,time(1):time(2)),ISClabels(time(1):time(2),:),ISCstruct.posIndices,tempIHClabels(time(1):time(2),:),IHCstruct.IHCcenters(newOrder,:),[],[0 max(img,[],'all')],h)
            elseif x == 1
                satisfied = 1;
            elseif x > 1
                    pkBinary = zeros(time(2)-time(1),size(activeIHCs,1));
                    %detect peaks in signal to create points for clustering
                    for j = 1:size(activeIHCs,1)
                     [pks,locs] = findpeaks(IHCrois(time(1):time(2),activeIHCs(j)),'MinPeakProminence',0.3);
                     pkBinary(locs,j) = 1;
                    end
                    [r,c] = find(pkBinary);
                    label = kmeans([r*20 c],x,'MaxIter',200); %time expanded so that this is the major clustering variable, rather than space
                    pkBinaryLab = pkBinary;
                    for j = 1:size(r,1)
                        pkBinaryLab(r(j),c(j)) = label(j);
                    end
                    figure; imagesc(pkBinaryLab);
                    pause(2);
                    %have to change raw and group labels, but only of second
                    %group
                    labelTimes = [];
                    for j = 1:x
                        [r,c] = find(pkBinaryLab == j);
                        labelTimes(j,:) = [j min(r) max(r)];
                    end
                    labelTimes = sortrows(labelTimes,2);

                    tempIHCgrouplabels = IHCgrouplabels; %to preserve original structure
                    tempIHCrawlabels = rawIHCgrouplabels;
                    labelsForEvent = [i];
                    for j = 2:x
                        tempBinary = pkBinaryLab == labelTimes(j,1);
                        tempBinary = imgaussfilt(single(tempBinary),[0.5 4]) > 0.1;
                        imagesc(tempBinary);
                        pause(2);
                        timeToChange = time(1) + labelTimes(j,2) - 2;
                        [r3,c3] = find(tempBinary); %smear across IHCs to get a more natural grouping
                        groupIHCs = rmoutliers(unique([c3]));
                        mask = tempIHCgrouplabels(timeToChange:time(2),activeIHCs(groupIHCs)) == i;
                        IHCgrouplabels(timeToChange:time(2),activeIHCs(groupIHCs)) = mask * (numLabels + 1);
                        rawIHCgrouplabels(timeToChange:time(2),activeIHCs(groupIHCs)) = mask * (max(rawIHCgrouplabels,[],'all') + 1);
                        numLabels = numLabels + 1;
                        labelsForEvent = [labelsForEvent; numLabels];
                    end

                    tempIHClabels = zeros(size(tempIHClabels));
                    for j = 1:size(labelsForEvent,1)
                        tempIHClabels = tempIHClabels + j*(IHCgrouplabels == labelsForEvent(j));
                    end
                    generateISCIHCActivityMovie(img(:,:,time(1):time(2)),ISClabels(time(1):time(2),:),ISCstruct.posIndices,tempIHClabels(time(1):time(2),:),IHCstruct.IHCcenters(newOrder,:),[],[0 max(img,[],'all')],h); 
                    satisfied = input('Are you satisfied (1) or not (0)?');
                    if ~satisfied %reset original values
                        IHCgrouplabels = tempIHCgrouplabels;
                        rawIHCgrouplabels = tempIHCrawlabels;
                    end
            end
        end
    end
end
%%

for i=1:max(IHCgrouplabels,[],'all')
    [r,c] = find(IHCgrouplabels == i);
    timeStart = min(r);
    timeEnd = max(r);
    tempISClabels = ISClabels(timeStart:timeEnd,:);
    tempISClabels = unique(tempISClabels);
    tempISClabels = tempISClabels(tempISClabels ~= 0);
    
    timeEnd = timeEnd + 5; %for correlations
    timeEnd(timeEnd > size(IHCrois,1)) = size(IHCrois,1);
    meanIHCsignal = mean(IHCrois(timeStart:timeEnd,c),2);
    
    if ~isempty(tempISClabels) & (timeStart ~= timeEnd)
        [r,c] = find(IHCgrouplabels == i);
        activeRange = [min(c) max(c)];
        %centroid measurements for overlapping IHC/ISC events
        centroidIHCs = mean(IHClocs(activeRange(1):activeRange(2),:),1);
        centroidISCs = [];
        meanISCsignal = [];
        for j=1:size(tempISClabels,1)
            [r,c] = find(ISClabels == tempISClabels(j));
            meanISCsignal(:,j) = mean(ISCrois(timeStart:timeEnd,c),2);
            
            ISCcs = unique(c);
            centroidISCs(j,:) = mean(ISCstruct.posIndices(ISCcs,1:2),1);
        end
        [b,idx] = min(pdist2(centroidIHCs,centroidISCs));
        centDist = pdist2(centroidIHCs,centroidISCs)';
        
        corrMat = corr([meanIHCsignal meanISCsignal]);
        corrMat = corrMat(2:end,1);
        comp = [[1:size(corrMat,1)]' centDist corrMat;];
        comp = comp(comp(:,2) < 100 & comp(:,3) > 0.2,:);
        [~,idx] = max(comp(:,3));
        idx = comp(idx,1);
        
        if ~isempty(idx)
            temp = rawIHCgrouplabels(IHCgrouplabels == i);
            temp = temp(temp~=0);
            temp = unique(temp);
            if ~isempty(temp)
                for j = 1:size(temp,1)
                    IHClabels(rawIHCgrouplabels == temp(j)) = tempISClabels(idx);
                end
            end
           % IHClabels(IHCgrouplabels == i) = tempISClabels(idx);
        else
            disp(['rejected, corr: ' num2str(corrMat(idx+1,1)) ' dist: ' num2str(b)]);
        end
   
    end
end
%%
 generateISCIHCActivityMovie(img,ISClabels,ISCstruct.posIndices,IHClabels,IHCstruct.IHCcenters(newOrder,:),[fp '\' name '_coordISCIHCaVid.mp4'],[0 max(img,[],'all')])

%%
%save the coordination
CoorStruct = struct();
CoorStruct.ISClabels = ISClabels;
CoorStruct.posIndices = ISCstruct.posIndices;
CoorStruct.IHClabels = IHClabels;
CoorStruct.IHCcenters = IHCstruct.IHCcenters;
save([fp '\' name '_coorStruct.mat'],'CoorStruct');