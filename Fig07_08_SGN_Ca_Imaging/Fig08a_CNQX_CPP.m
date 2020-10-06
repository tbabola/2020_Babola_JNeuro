%% P0 base WT CNQX/CPP Experiments
clear all; close all;

dir = './Data/P0 WT*/base/*struct2*.mat';
files = findPairs(dir,'cpp');
SGNstructs = loadCellStructs(files);

for i=1:2:size(SGNstructs,2)
    if size([SGNstructs(i).events.isEmpty]) == size([SGNstructs(i+1).events.isEmpty])
        masks = [SGNstructs(i).events.isEmpty] + [SGNstructs(i+1).events.isEmpty];
        masks = masks >= 2;
        SGNstructs(i).meanFreq_ActiveD = mean([SGNstructs(i).events(~masks).frequency]);
        SGNstructs(i+1).meanFreq_ActiveD = mean([SGNstructs(i+1).events(~masks).frequency]);
    else
        i
    end
end

ROIsForGroup = 35;
GroupBase = groupedActivity(SGNstructs,ROIsForGroup);
% for i = 1:50
%     GroupBase = groupedActivity(SGNstructs,i);
%     freqsWT(:,i) = [GroupBase(1:2:end).freq];
%     freqsCPP(:,i) = [GroupBase(2:2:end).freq];
%     roisWT(:,i) = [GroupBase(1:2:end).meanROIs];
%     roisCPP(:,i) = [GroupBase(2:2:end).meanROIs];
% end
%
% figure; plot(mean(freqsWT,1)); hold on;
% plot(mean(freqsCPP,1));
% 
% figure; plot(mean(roisWT,1)); hold on;
% plot(mean(roisCPP,1));

conSGNstructs = SGNstructs(1:2:end);
cppSGNstructs = SGNstructs(2:2:end);
CNQXcolor = [102 153 255]/255;
h = compare2P([conSGNstructs.meanCorr_Active]', [cppSGNstructs.meanCorr_Active]', {"Baseline","+CNQX/CPP"}, 'Correlation coefficient', [3 2], [15 10],'k',CNQXcolor);
xtickangle(45);
ylim([0 1]);
yticks(0:0.25:1);
h1 =compare2P([conSGNstructs.meanFreq_ActiveD]', [cppSGNstructs.meanFreq_ActiveD]', {"Baseline","+CNQX/CPP"}, 'Active ROI frequency', [3 2], [15 10],'k',CNQXcolor);
xtickangle(45);
ylim([0 0.5]);
yticks(0:0.1:.5);
h2 = compare2P([conSGNstructs.meanAmplitude_Active]', [cppSGNstructs.meanAmplitude_Active]', {"Baseline","+CNQX/CPP"}, 'Amp', [3 2], [15 10],'k',CNQXcolor);
xtickangle(45);
ylim([0 3]);
yticks(0:1:3);
h3 = compare2P([conSGNstructs.activeArea]', [cppSGNstructs.activeArea]', {"Baseline","+CNQX/CPP"}, 'Active area (%)', [3 2], [15 10],'k',CNQXcolor);
xtickangle(45);
ylim([0 1]);
yticks(0:.25:1);
h4 = compare2P([conSGNstructs.meanHW_Active]', [cppSGNstructs.meanHW_Active]', {"Baseline","+CNQX/CPP"}, 'half-width (s)', [3 2], [15 10],'k',CNQXcolor);
xtickangle(45);
ylim([0 10]);
yticks(0:2.5:10);
h5 = compare2P([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]', {"Baseline","+CNQX/CPP"}, 'Correlated events per min', [3 2], [15 10],'k',CNQXcolor);
xtickangle(45);
ylim([0 3]);
yticks(0:1:3);
%handleTheSubplot({h5,h1,h2,h3,h4,h},[1 6]);
handleTheSubplot({h5,h,h3,h1},[1 4])
figQuality(gcf,gca,[4.2 1.75]);
export_fig('.\EPS Panels\P0_WT_CNQX_CPP.eps');



%%

disp(['n = ' num2str(size(SGNstructs,2)/2)]);
[~,p, ci, stats] = ttest([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]');
disp(['Correlated events per minute p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

[~,p,~,stats] = ttest([conSGNstructs.meanCorr_Active]', [cppSGNstructs.meanCorr_Active]');
disp(['Corr p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

[~,p,~,stats] = ttest([conSGNstructs.activeArea]', [cppSGNstructs.activeArea]');
disp(['Active Roi area (%) p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

[~,p,~,stats] = ttest([conSGNstructs.meanFreq_ActiveD]', [cppSGNstructs.meanFreq_ActiveD]');
disp(['Freq p-value = ' num2str(p) '  t-stat: ' num2str(stats.tstat) '    tdf: ' num2str(stats.df)]);

%not used
% [~,p] = ttest([conSGNstructs.meanAmplitude_Active]', [cppSGNstructs.meanAmplitude_Active]');
% disp(['Amplitude p-value = ' num2str(p)]);
% 
% [~,p] = ttest([conSGNstructs.meanHW_Active]', [cppSGNstructs.meanHW_Active]');
% disp(['HW p-value = ' num2str(p)]);


numTests = 6;
disp(['Bonferroni correction = 0.05 / ' num2str(numTests) ' = ' num2str(0.05/numTests)]);

%%
%cochlea_5_base_2 and base_3_cnqx.cpp_2
%baseline 601-875, MRS 301-575
%median filt 3d (1x1x1) then normalize contrast (18 to 3400), in
%data/Time_color_Processed

cellNum = 4
rois = conSGNstructs(cellNum).rois;
roisCNQX = cppSGNstructs(cellNum).rois;
[t,n] = size(rois);

numToShow = 100;
rng(1)
temp = randperm(n)';
temp = sort(temp(1:numToShow));
isEmptyCon = [conSGNstructs(cellNum).events(temp).isEmpty];
isEmptyCpp = [cppSGNstructs(cellNum).events(temp).isEmpty];
offsets = -0.2*repmat(1:numToShow,t,1);

for i = 1:numToShow
    rois(:,i) = smooth(rois(:,temp(i)),3);
    roisCNQX(:,i) = smooth(roisCNQX(:,temp(i)),3);
end

figure; plot(rois(1:275,isEmptyCon)+offsets(1:275,isEmptyCon),'Color',[0.7 0.7 0.7]); hold on;
plot(rois(1:275,~isEmptyCon)+offsets(1:275,~isEmptyCon),'Color','k');
ylim([-21 5]);
xlim([0 275]);
figQuality(gcf,gca,[2.4 1.6])
export_fig('.\EPS Panels\baseline_CNQX_example.eps')
% 
figure; plot(roisCNQX(:,isEmptyCpp)+offsets(1:275,isEmptyCpp),'Color',[0.7 0.7 0.7]); hold on;
plot(roisCNQX(:,~isEmptyCpp)+offsets(1:275,~isEmptyCpp),'Color','k');
ylim([-21 5]);
xlim([0 275]);
figQuality(gcf,gca,[2.4 1.6])
export_fig('.\EPS Panels\CNQX_example.eps')

%%
%% E165 base WT CNQX/CPP Experiments
clear all; close all;

dir = './Data/E16.5*/base/*struct2*.mat';
files = findPairs(dir,'cpp');
SGNstructs = loadCellStructs(files);

for i=1:2:size(SGNstructs,2)
    if size([SGNstructs(i).events.isEmpty]) == size([SGNstructs(i+1).events.isEmpty])
        masks = [SGNstructs(i).events.isEmpty] + [SGNstructs(i+1).events.isEmpty];
        masks = masks >= 2;
        SGNstructs(i).meanFreq_ActiveD = mean([SGNstructs(i).events(~masks).frequency]);
        SGNstructs(i+1).meanFreq_ActiveD = mean([SGNstructs(i+1).events(~masks).frequency]);
    else
        i
    end
end

ROIsForGroup = 35;
GroupBase = groupedActivity(SGNstructs,ROIsForGroup);
% for i = 1:50
%     GroupBase = groupedActivity(SGNstructs,i);
%     freqsWT(:,i) = [GroupBase(1:2:end).freq];
%     freqsCPP(:,i) = [GroupBase(2:2:end).freq];
%     roisWT(:,i) = [GroupBase(1:2:end).meanROIs];
%     roisCPP(:,i) = [GroupBase(2:2:end).meanROIs];
% end
%
% figure; plot(mean(freqsWT,1)); hold on;
% plot(mean(freqsCPP,1));
% 
% figure; plot(mean(roisWT,1)); hold on;
% plot(mean(roisCPP,1));

conSGNstructs = SGNstructs(1:2:end);
cppSGNstructs = SGNstructs(2:2:end);

h = compare2P([conSGNstructs.meanCorr_Active]', [cppSGNstructs.meanCorr_Active]', {"WT Baseline","+CNQX/CPP"}, 'Corr', [3 2], [15 10]);
ylim([0 1]);
yticks(0:0.25:1);
h1 =compare2P([conSGNstructs.meanFreq_ActiveD]', [cppSGNstructs.meanFreq_ActiveD]', {"WT Baseline","+CNQX/CPP"}, 'freq', [3 2], [15 10]);
ylim([0 0.5]);
yticks(0:0.1:.5);
h2 = compare2P([conSGNstructs.meanAmplitude_Active]', [cppSGNstructs.meanAmplitude_Active]', {"WT Baseline","+CNQX/CPP"}, 'Amp', [3 2], [15 10]);
ylim([0 3]);
yticks(0:1:3);
h3 = compare2P([conSGNstructs.activeArea]', [cppSGNstructs.activeArea]', {"WT Baseline","+CNQX/CPP"}, 'Area', [3 2], [15 10]);
ylim([0 1]);
yticks(0:.25:1);
h4 = compare2P([conSGNstructs.meanHW_Active]', [cppSGNstructs.meanHW_Active]', {"WT Baseline","+CNQX/CPP"}, 'half-width (s)', [3 2], [15 10]);
ylim([0 10]);
yticks(0:2.5:10);
h5 = compare2P([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]', {"WT Baseline","+CNQX/CPP"}, 'Correlated events per minute', [3 2], [15 10]);
ylim([0 3]);
yticks(0:1:3);
%handleTheSubplot({h5,h1,h2,h3,h4,h},[1 6]);
handleTheSubplot({h5,h,h3,h1},[1 4])
figQuality(gcf,gca,[4 1.25]);
export_fig('.\EPS Panels\E165_WT_CNQX_CPP.eps');

[~,p] = ttest([conSGNstructs.meanCorr_Active]', [cppSGNstructs.meanCorr_Active]');

disp(['n = ' num2str(size(SGNstructs,2)/2)]);
disp(['Corr p-value = ' num2str(p)]);
[~,p] = ttest([conSGNstructs.meanFreq_ActiveD]', [cppSGNstructs.meanFreq_ActiveD]');
disp(['Freq p-value = ' num2str(p)]);
[~,p] = ttest([conSGNstructs.meanAmplitude_Active]', [cppSGNstructs.meanAmplitude_Active]');
disp(['Amplitude p-value = ' num2str(p)]);
[~,p] = ttest([conSGNstructs.activeArea]', [cppSGNstructs.activeArea]');
disp(['Active Roi area (%) p-value = ' num2str(p)]);
[~,p] = ttest([conSGNstructs.meanHW_Active]', [cppSGNstructs.meanHW_Active]');
disp(['HW p-value = ' num2str(p)]);
[~,p] = ttest([GroupBase(1:2:end).freq]', [GroupBase(2:2:end).freq]');
disp(['Correlated event p-value = ' num2str(p)]);

numTests = 6;
disp(['Bonferroni correction = 0.05 / ' num2str(numTests) ' = ' num2str(0.05/numTests)]);

%%
cellNum = 2;
rois = conSGNstructs(cellNum).rois;
roisCNQX = cppSGNstructs(cellNum).rois;
[t,n] = size(rois);

numToShow = 100;
rng(1)
temp = randperm(n)'
temp = sort(temp(1:numToShow));
for i = 1:numToShow
    rois(:,i) = smooth(rois(:,i),3);
    roisCNQX(:,i) = smooth(roisCNQX(:,i),3);
end

figure; plot(rois(1:300,temp) - 0.2*repmat(1:numToShow,300,1),'Color','k');
ylim([-20 5]);
figQuality(gcf,gca,[2.6 1.5])
export_fig('.\EPS Panels\E165_baseline_CNQX_example.eps')

figure; plot(roisCNQX(:,temp) - 0.2*repmat(1:numToShow,size(roisCNQX,1),1),'Color','k');
ylim([-20 5])
figQuality(gcf,gca,[2.6 1.5])
export_fig('.\EPS Panels\E165_CNQX_example.eps')