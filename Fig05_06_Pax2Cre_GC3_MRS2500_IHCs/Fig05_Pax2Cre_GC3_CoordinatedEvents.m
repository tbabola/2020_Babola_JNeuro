%% Figure 5: Pax2-Cre;R26-lsl-GCaMP3 + MRS2500
%% 
% Panel A-D is image 45, std projection of frames 775-850 made in imageJ, raw file is ./Data/P0 MRS2500/PanelA.tif
% and straightened in ImageJ using roi ./Data/P0 MRS2500/Image45_staightenROI.roi.
% Panels C/D was made by using median filter over Image45 (straightened;
% 1x1x1) and using Time-Color plot in ImageJ.
%% Panel E (traces)
%ISCs
ISC = load(ISCfiles{5});
ISCstruct = ISC.ISCstruct;
rois = ISCstruct.rois;
numToShow = 75; n = size(rois,2)
if numToShow > n
    numToShow = n;
end
temp = randperm(n)';
temp = sort(temp(1:numToShow));
figure; plotOffset(rois(1:1800,temp),0.5)
figQuality(gcf,gca,[4.25 2]);
export_fig('.\EPS Panels\PanelE_ISCs.eps');

%IHCs
IHCfiles = loadFileList('.\Data\P0 MRS2500\*_IHCstructwhole*');
ISCfiles = loadFileList('.\Data\P0 MRS2500\*_ISC*');
IHC = load(IHCfiles{1})
IHCstruct = IHC.IHCstruct;

rois = IHCstruct.smoothroi;
numToShow = 75; n = size(rois,2)
if numToShow > n
    numToShow = n;
end
temp = randperm(n)';
temp = sort(temp(1:numToShow));

figure; plotOffset(IHCstruct.smoothroi(1:1800,temp),0.5)
figQuality(gcf,gca,[4.25 2]);
export_fig('.\EPS Panels\PanelE_IHCs.eps');
%% Panel F
%Quantification
IHCstructs = loadCellStructs('.\Data\P0 MRS2500\*_IHCstruct.mat');
IHCstructsMRS = loadCellStructs('.\Data\P0 MRS2500\*_IHCstructMRS2500.mat');
ISCstructs = loadCellStructs(['.\Data\P0 MRS2500\*_ISCdata.mat*']);

IHCfreq = []; IHCfreqMRS = [];
ISCfreq = []; ISCfreqMRS = [];

for i = 1:size(IHCstructs,2)
    IHCfreq(i) = size(IHCstructs(i).event,2)/5;
    IHCfreqMRS(i) = size(IHCstructsMRS(i).event,2)/5;
    IHCarea(i) = mean([IHCstructs(i).event.area]);
    if IHCfreqMRS(i) ~= 0
        [IHCstructsMRS(i).event.area]
        IHCareaMRS(i) = mean([IHCstructsMRS(i).event.area]);
    else
        IHCareaMRS(i) = 0;
    end
    IHCcorr(i) = mean(prctile(IHCstructs(i).corr,90));
    IHCcorrMRS(i) = mean(prctile(IHCstructsMRS(i).corr,90));
    
    tempISC = [ISCstructs(i).event([ISCstructs(i).event.timeStart] <= 600)];
    tempISCMRS = [ISCstructs(i).event([ISCstructs(i).event.timeStart] > 1200 & [ISCstructs(i).event.timeStart] <=1800)];
    ISCfreq(i) = size(tempISC,2)/5;
    ISCfreqMRS(i) = size(tempISCMRS,2)/5;
    ISCarea(i) = mean([tempISC.area]);
    if ISCfreqMRS(i) ~= 0
        ISCareaMRS(i) = mean([tempISCMRS.area]);
    else
        ISCareaMRS(i) = 0;
    end
end

%% Panel F Plots
h1 = compare2Px2(IHCfreq', IHCfreqMRS', ISCfreq', ISCfreqMRS',{'Baseline','MRS2500','Baseline','MRS2500'}, 'Events per minute', [1.5 2], [15 8], 'k', 'r')
xtickangle(45);
ylim([0 15]);
figQuality(gcf,gca,[1.5 2]);

h2 = compare2Px2(IHCarea', IHCareaMRS', ISCarea', ISCareaMRS',{'Baseline','MRS2500','Baseline','MRS2500'}, '# ROIs per event', [1.5 2], [15 8], 'k', 'r')
xtickangle(45);
ylim([0 15]);
figQuality(gcf,gca,[1.5 2]);

h3 = compare2P(IHCcorr',IHCcorrMRS',{'Baseline','MRS2500'},'Correlation Coef',[3 2],[15 10],'k','r');
xtickangle(45);
ylim([0 1]);
figQuality(gcf,gca,[1.5 2]);

handleTheSubplot({h1,h2,h3},[1 3]);
figQuality(gcf,gca,[3.9,1.25]);
export_fig('./EPS Panels/PanelE.eps')

%% Statistics
[h,p,ci,stats] = ttest(IHCfreq,IHCfreqMRS)
[h,p2,ci2,stats2] = ttest(ISCfreq,ISCfreqMRS)

[h,p3,ci3,stats3] = ttest(IHCarea,IHCareaMRS)
[h,p4,ci4,stats4] = ttest(ISCarea,ISCareaMRS)

[h,p5,ci5,stats5] = ttest(IHCcorr,IHCcorrMRS)