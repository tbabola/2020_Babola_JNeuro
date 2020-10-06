%% Figure 6: IHC events in Pax2-Cre;R26-lsl-GCaMP3 mice

% Panel A is image 45, std projection of frames 318-322 made in imageJ, raw file is ./Data/P0 MRS2500/PanelA.tif
% and straightened in ImageJ using roi ./Data/P0 MRS2500/Image45_staightenROI.roi.
% Panels B/C was made by using median filter over Image45 (straightened;
% 1x1x1) and using Time-Color plot in ImageJ from 0 to 600 frames (5 minutes).

%Panel A (middle) is image 41.
%Panel A (base) image 39.
%% load the base data
addpath(genpath('..\MATLAB Functions'));
ISCstructs = loadCellStructs('.\Data\P0 Base\*ISCdata.mat');
IHCstructs = loadCellStructs('.\Data\P0 Base\*IHCstruct.mat');

structs = loadCellStructs('.\Data\P0 Base\*coorStruct2.mat');
structs = getAreas(structs);
structs = getAreasAll(structs,ISCstructs);
BaseSmallEvents = getSmallEvents(structs);
%% Panel B, base
[~,slopes] = plotScatter(structs);
disp([num2str(mean(slopes)) ' +/- ' num2str(sterr(slopes,2)) ' IHC/ISC' ]);
ylim([0 80]); xlim([0 80]);
xlabel('# of ISC ROIs activated');
ylabel('# of IHCs activated');
figQuality(gcf,gca,[1.5 1.3]);
export_fig('.\EPS Panels\BaseEvents.eps');

%some stats
[means] = statScatter(structs);
m = mean(means,1);
s = sterr(means,1);
disp(['Base ISCs active: ' num2str(m(1)) ' +/- ' num2str(s(1))]);
disp(['Base IHCs active: ' num2str(m(2)) ' +/- ' num2str(s(2))]);
base_m = means;
base_s = s;
%% Panel C, base
IHCsignals = plotCenteredSignals(structs,IHCstructs);
xlim([0 800]);
ylim([0 2.5]);
figQuality(gcf,gca,[3 1]);
print('.\EPS Panels\BaseCentered.pdf','-dpdf');
disp([num2str(size(IHCsignals,2)) ' events.']);
%get mean + sem of middle cell
[means,sem] = grpstats(max(IHCsignals(:,:,21)),ones(size(IHCsignals,2),1),{'mean','sem'})
IHCb = IHCsignals;
%% load the apex data
IHCstructs = loadCellStructs('.\Data\P0 Apex Apical\*IHCstruct.mat')
ISCstructs = loadCellStructs('.\Data\P0 Apex Apical\*ISCdata.mat')
structs = loadCellStructs('.\Data\P0 Apex Apical\*coorStruct2.mat');
structs = getAreas(structs);
structs = getAreasAll(structs,ISCstructs);
ApexsmallEvents = getSmallEvents(structs);
%% Panel B, apex
[~,slopes] = plotScatter(structs);
disp([num2str(mean(slopes)) ' +/- ' num2str(sterr(slopes,2)) ' IHC/ISC' ]);
ylim([0 80]); xlim([0 80]);
xlabel('# of ISC ROIs activated');
ylabel('# of IHCs activated');
figQuality(gcf,gca,[1.5 1.3]);
export_fig('.\EPS Panels\ApexEvents.eps');

%some stats
[means] = statScatter(structs);
m = mean(means,1);
s = sterr(means,1);
disp(['Apex ISCs active: ' num2str(m(1)) ' +/- ' num2str(s(1))]);
disp(['Apex IHCs active: ' num2str(m(2)) ' +/- ' num2str(s(2))]);
apex_m = means;
apex_s = s;
%% Panel C, apex
IHCsignals = plotCenteredSignals(structs,IHCstructs);
xlim([0 800]);
ylim([0 2.5]);
figQuality(gcf,gca,[3 1]);
print('.\EPS Panels\ApexCentered.pdf','-dpdf');
disp([num2str(size(IHCsignals,2)) ' events.']);
%get mean + sem of middle cell
[means,sem] = grpstats(max(IHCsignals(:,:,21)),ones(size(IHCsignals,2),1),{'mean','sem'})
IHCa = IHCsignals;
%% load the middle data
IHCstructs = loadCellStructs('.\Data\P0 Apex Middle\*IHCstruct.mat')
ISCstructs = loadCellStructs('.\Data\P0 Apex Middle\*ISCdata.mat')
structs = loadCellStructs('.\Data\P0 Apex Middle\*coorStruct2.mat');
structs = getAreas(structs);
structs = getAreasAll(structs,ISCstructs);
MiddleSmallEvents = getSmallEvents(structs);
%% Panel B, middle
[~,slopes] = plotScatter(structs);
disp([num2str(mean(slopes)) ' +/- ' num2str(sterr(slopes,2)) ' IHC/ISC' ]);
ylim([0 80]); xlim([0 80]);
xlabel('# of ISC ROIs activated');
ylabel('# of IHCs activated');
figQuality(gcf,gca,[1.5 1.3]);
export_fig('.\EPS Panels\MiddleEvents.eps');
%some stats
[means] = statScatter(structs);
m = mean(means,1);
s = sterr(means,1);
disp(['Middle ISCs active: ' num2str(m(1)) ' +/- ' num2str(s(1))]);
disp(['Middle IHCs active: ' num2str(m(2)) ' +/- ' num2str(s(2))]);
m_m = means;
m_s = s;
%% Panel C, middle
IHCsignals = plotCenteredSignals(structs,IHCstructs);
xlim([0 800]);
ylim([0 2.5]);
figQuality(gcf,gca,[3 1]);
print('.\EPS Panels\MiddleCentered.pdf','-dpdf');
disp([num2str(size(IHCsignals,2)) ' events.']);
%get mean + sem of middle cell
[means,sem] = grpstats(max(IHCsignals(:,:,21)),ones(size(IHCsignals,2),1),{'mean','sem'})
IHCm = IHCsignals;
%% Amplitudes
ISCstructs = loadCellStructs('.\Data\P0 Base\newISCstructs\*ISCdata.mat');
[base_amps, mb_amps, base_mamps, mb_mamps] = getAmps(ISCstructs);
ISCstructs = loadCellStructs('.\Data\P0 Apex Middle\newISCstructs\*ISCdata.mat');
[mid_amps, mm_amps, m_mamps, mm_mamps] = getAmps(ISCstructs);
ISCstructs = loadCellStructs('.\Data\P0 Apex Apical\newISCstructs\*ISCdata.mat');
[apex_amps, ma_amps, a_mamps, ma_mamps] = getAmps(ISCstructs);

compare3(mb_amps,mm_amps,ma_amps,{'Base','Middle','Apex'}, "Max Amplitude", [3 2])
compare3(mb_mamps,mm_mamps,ma_mamps,{'Base','Middle','Apex'}, "Mean Amplitude", [3 2])

anova1([mb_amps,mm_amps,ma_amps]',[ones(size(mb_amps,2),1); 2*ones(size(mm_amps,2),1); 3*ones(size(ma_amps,2),1)])
anova1([mb_mamps,mm_mamps,ma_mamps]',[ones(size(mb_mamps,2),1); 2*ones(size(mm_mamps,2),1); 3*ones(size(ma_mamps,2),1)]   )
[means, std] = grpstats([mb_mamps,mm_mamps,ma_mamps]',[ones(size(mb_mamps,2),1); 2*ones(size(mm_mamps,2),1); 3*ones(size(ma_mamps,2),1)],{'mean','sem'});
%%
%%Stats on event size
tempAmps = [max(IHCa(:,:,21)) max(IHCm(:,:,21)) max(IHCb(:,:,21))]';
tempGroups =  [ones(size(max(IHCa(:,:,21)),2),1); 2*ones(size(max(IHCm(:,:,21)),2),1); 3*ones(size(max(IHCb(:,:,21)),2),1)];
[p,~,stats] = anova1(tempAmps,tempGroups);
c = multcompare(stats);

%stats on ISCs activated
tempISCs = [apex_m(:,1); m_m(:,1); base_m(:,1)];
tempGroups = [ones(size(apex_m(:,1),1),1); 2*ones(size(m_m(:,1),1),1); 3*ones(size(base_m(:,1),1),1)];
p = anova1(tempISCs, tempGroups)



%%
function [maxAmps, meanMaxAmps, meanAmps, meanMeanAmps] = getAmps(ISCstructs)
    maxAmps = [];
    meanMaxAmps = [];
    meanAmps = [];
    meanMeanAmps = [];
    for i = 1:size(ISCstructs,2)
        temp = [ISCstructs(i).event.maxAmplitude]';
        maxAmps = [maxAmps; temp];
        meanMaxAmps(i) = mean(temp,1);
        
        temp = [ISCstructs(i).event.meanAmplitude]';
        meanAmps = [meanAmps ; temp];
        meanMeanAmps(i) = mean(temp,1);
        
        
    end
end

function structs = getAreas(structs)
    for i=1:size(structs,2)
        ISClabels = structs(i).ISClabels;
        IHClabels = structs(i).IHClabels;

        events = struct();
        for j = 1:max(ISClabels,[],'all')
           [~,c] = find(ISClabels == j);
           events(j).areaISC = size(unique(c),1);
           [~,c] = find(IHClabels == j);
           events(j).areaIHC = size(unique(c),1);
           if ~isempty(unique(c))
               events(j).areaIHC = size(unique(c),1);
           else
               events(j).areaIHC = 0;
           end
        end
        structs(i).events = events;
    end
end

function structs = getAreasAll(structs,ISCstructs)
    for i=1:size(structs,2)
        multiISClabels = ISCstructs(i).labelRoi;
        ISClabels = structs(i).ISClabels;
        IHClabels = structs(i).IHClabels;

        eventsAll = struct();
        for j = 1:max([ISClabels IHClabels],[],'all')
           [r,c] = find(ISClabels == j);
           
           if ~isempty(unique(c))
             eventsAll(j).areaISC = size(unique(c),1);
             multilabel = unique(multiISClabels(ISClabels == j));
             multilabel = multilabel(multilabel ~=0);
             if ~isempty(multilabel) & size(multilabel,1) > 1
                 eventsAll(j).multiISCs = size(multilabel,1);
                 eventsAll(j).maxResponse = 0; 
                 eventsAll(j).avgResponse = 0;
             else
                 eventsAll(j).multiISCs = 0;
                 temp = ISCstructs(i).rois;
                 temp = temp(min(r):max(r),c);
                 eventsAll(j).maxResponse = max(temp,[],'all');
                 eventsAll(j).avgResponse = max(mean(temp,2),[],'all');
             end
           else
               eventsAll(j).areaISC = 0;
               eventsAll(j).multiISCs = 0;
               eventsAll(j).maxResponse = 0;  
               eventsAll(j).avgResponse = 0;
           end
           
           
           [~,c] = find(IHClabels == j);
           if ~isempty(unique(c))
               eventsAll(j).areaIHC = size(unique(c),1);
           else
               eventsAll(j).areaIHC = 0;
           end
        end
        structs(i).eventsAll = eventsAll;
    end
end

function [h, slopes] = plotScatter(structs)
    h = figure; hold off;
    for i=1:size(structs,2)
        mask = [structs(i).eventsAll.multiISCs] == 0 & [structs(i).eventsAll.areaISC] > 0;
        tempx = [structs(i).eventsAll.areaISC];
        tempy = [structs(i).eventsAll.areaIHC];
        tempy = tempy(mask);
        tempx = tempx(mask);

        scatter(tempx',tempy',20,[0.7 0.7 0.7],'.','jitter','on');
        f = fit(tempx',tempy','a*x');
        slopes(i) = f.a;
        x = 0:80;
        plot(x,f(x),'Color',[0.7 0.7 0.7]);
        if i == 1
            hold on;
        end
    end
    
     for i = 1:size(structs,2)
        mask = [structs(i).eventsAll.multiISCs] == 0 & [structs(i).eventsAll.areaISC] > 0;
        tempx = [structs(i).eventsAll.areaISC];
        tempy = [structs(i).eventsAll.areaIHC];
        tempy = tempy(mask);
        tempx = tempx(mask);
        yerror = sterr(tempy,2);
        xerror = sterr(tempx,2);
        errorbar(mean(tempx,2),mean(tempy,2),yerror,yerror,... 
            xerror, xerror,'CapSize',0,'Color','k','Marker','.','MarkerSize',7);
    end
end

function means = statScatter(structs)
    for i=1:size(structs,2)
        mask = [structs(i).eventsAll.multiISCs] == 0 & [structs(i).eventsAll.areaISC] > 0;
        tempx = [structs(i).eventsAll.areaISC];
        tempy = [structs(i).eventsAll.areaIHC];
        tempy = tempy(mask);
        tempx = tempx(mask);
        
        means(i,:) = mean([tempx' tempy'],1);
    end
end

function smallEvents = getSmallEvents(structs)
    smallEvents = [];
    for i=1:size(structs,2)
        ISCarea = [structs(i).eventsAll.areaISC];
        IHCarea = [structs(i).eventsAll.areaIHC];
        smallEventCount = 0;
        for j = 1:size(IHCarea,2)
            if ISCarea(j) == 0
                [~,c] = find(structs(i).IHClabels == j);
                IHCc = unique(c);
                IHCc = IHCc(IHCc~=0);
                if size(IHCc,1) < 2
                    smallEventCount = smallEventCount + 1;
                end
            end
        end
        smallEvents(i) = smallEventCount;
    end
end

function locs = getIHClocations(structs)
    locs = [];
    for i = 1:size(structs,2)
        tempISC = [structs(i).eventsAll.areaISC];
        tempIHC = [structs(i).eventsAll.areaIHC];
        for j = 1:size(tempISC,2)
            [r2,c2] = find([structs(i).IHClabels structs(i).ISClabels] == j);
            if ~isempty(r2)
                time = min(r2);
            else
                time = NaN;
            end
            [r,c] = find(structs(i).IHClabels == j);
            min(r2);
            c = c/size(structs(i).IHClabels,2);
             if tempISC(j) == 0
                if ~isempty(c)
                    locs = [locs; min(c) max(c) 1 i j time];
                else
                   locs = [locs; 0 0 1 i j time]; 
                end
            else
                if ~isempty(c)
                    locs = [locs; min(c) max(c) 0 i j time];
                else
                    locs = [locs; 0 0 0 i j time]; 
                end
            end
        end
    end
end

function IHCsignals = plotCenteredSignals(structs, IHCstructs)
    numSide = 20;
    IHCsignals = zeros(0,numSide*2+1); counts = zeros(1,numSide*2+1);
    totalCount = 1;
    for i=1:size(structs,2)
        ISCa = [structs(i).eventsAll.areaISC];
        IHCa = [structs(i).eventsAll.areaIHC];

        mask = IHCa > 10 & ISCa > 0;
        idx = find(mask);
        for j = 1:size(idx,2)
            [r,c] = find(structs(i).IHClabels == idx(j));
            [~,cISC] = find(structs(i).ISClabels == idx(j));
            cISC = unique(cISC);
            ISCcentroid = mean(structs(i).posIndices(cISC,1:2),1);
            c = [min(c)-1:max(c)+1]';
            c = c(c > 0 & c <= size(structs(i).IHClabels,2));
            IHCcenters = structs(i).IHCcenters(c,:);
            [~,minDist] = min(pdist2(ISCcentroid,IHCcenters));
            timeStart = min(r) - 2;
            timeStart(timeStart < 1) = 1;
            timeEnd = timeStart + 15;
            timeEnd(timeEnd > 600) = 600;

            IHCsignals([timeStart:timeEnd]-timeStart+1,totalCount,numSide+1) = IHCstructs(i).roisignalNorm(timeStart:timeEnd,c(minDist));

            for k=1:numSide
                if minDist-k > 1
                    IHCsignals([timeStart:timeEnd]-timeStart+1,totalCount,numSide+1-k) = IHCstructs(i).roisignalNorm(timeStart:timeEnd,c(minDist-k));
                else
                    IHCsignals([timeStart:timeEnd]-timeStart+1,totalCount,numSide+1-k) =0;
                end
                if minDist+k < size(c,1)
                    IHCsignals([timeStart:timeEnd]-timeStart+1,totalCount, numSide+1+k) = IHCstructs(i).roisignalNorm(timeStart:timeEnd,c(minDist+k));
                else
                    IHCsignals([timeStart:timeEnd]-timeStart+1,totalCount, numSide+1+k) = 0;
                end
            end
            totalCount = totalCount + 1;

        end
    end
    
    meanSig = squeeze(nanmean(IHCsignals,2));
    meanSig = [meanSig; NaN([3 size(meanSig,2)])];
    stdSig = squeeze(nanstd(IHCsignals,[],2))/sqrt(size(IHCsignals,2));
    stdSig = [stdSig; NaN([3 size(meanSig,2)])];
    [~,locs] = findpeaks(meanSig(:),'MinPeakDistance',10);
    x = 1:size(meanSig(:));
    xshade = x(locs);
    meanshade = meanSig(locs);
    errshade = stdSig(locs);
    figure;
    shadedErrorBar(xshade,meanshade,errshade); hold on;
    plot(meanSig(:),'k');
end