%% 1 animal 191126, 2 animals 191203, 3 animals 191204, 2 animals 191205
conditionTimes = [0 300; 480 780; 2100 2400];

crens = loadCellStructs('.\Data\E16\*.mat');
[freqs,area,aream] = getStats(crens,conditionTimes);
dataBlock = table(repmat(1,size(freqs,1),1),freqs,area,aream);

conditionTimes = [0 600; 900 1500; 2400 3000];
crens = loadCellStructs('.\Data\P0-2\*.mat');
[freqs,area,aream] = getStats(crens,conditionTimes);
dataBlock = [dataBlock; table(repmat(1.7,size(freqs,1),1),freqs,area,aream)];

crens7_8 = loadCellStructs('.\Data\P7-8\*.mat');
[freqs,area,aream] = getStats(crens7_8,conditionTimes);
dataBlock = [dataBlock; table(repmat(3,size(freqs,1),1),freqs,area,aream)];

crens10_12 = loadCellStructs('.\Data\P10-12\*.mat');
[freqs,area,aream] = getStats(crens10_12,conditionTimes);
dataBlock = [dataBlock; table(repmat(4,size(freqs,1),1),freqs,area,aream)];
dataStats = grpstats(dataBlock,'Var1',{'mean','sem'});

crens_0_base = loadCellStructs('.\Data\P0-2\base\*.mat');
[freqs,area,aream] = getStats(crens_0_base,conditionTimes);
%prep 4 washout did not happen because of air bubble and bath overflow
freqs(4,3) = NaN;
area(4,3) = NaN;
aream(4,3) = NaN;
dataBlock = [dataBlock; table(repmat(2.3,size(freqs,1),1),freqs,area,aream)];
dataStats = grpstats(dataBlock,'Var1',{'mean','sem'});

%%
meanMarkSize = 12;
color1 = 'k';
figure;
offset = 0.2;
points = [-offset 0 offset];
xs = repmat(dataBlock.Var1,1,3);
xs = xs+points;
ys = dataBlock.freqs;
plot(xs',ys','-','Color',[0.7 0.7 0.7]);
ylim([0 5]);
yticks(0:2.5:5);
hold on;
errorbar([1 1.7 2.3 3 4] - offset, dataStats.mean_freqs(:,1), dataStats.sem_freqs(:,1),'LineStyle', 'none','LineWidth',1,'Color',color1,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
errorbar([1 1.7 2.3 3 4], dataStats.mean_freqs(:,2), dataStats.sem_freqs(:,2),'LineStyle', 'none','LineWidth',1,'Color','r','CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
errorbar([1 1.7 2.3 3 4] + offset, dataStats.mean_freqs(:,3), dataStats.sem_freqs(:,3),'LineStyle', 'none','LineWidth',1,'Color',[0.4 0.4 0.4],'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
ylabel("Crenations per minute");
xticks(1:4);
xticklabels({'E16-17','P0-2','P6-8','P10-12'});

figQuality(gcf,gca,[3.75 1.25]);
export_fig('.\EPS Panels\crenation_Frequency.eps');
% meanMarkSize = 12;
% color1 = 'k';
% figure;
% offset = 0.3;
% points = [-offset 0 offset];
% ys = dataBlock.aream;
% plot(xs',ys','-','Color',[0.7 0.7 0.7]);
% ylim([0 5]);
% yticks(0:2.5:5);
% hold on;
% errorbar([1 2 3 4] - offset, dataStats.mean_area(:,1), dataStats.sem_area(:,1),'LineStyle', 'none','LineWidth',1,'Color',color1,'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
% errorbar([1 2 3 4], dataStats.mean_area(:,2), dataStats.sem_area(:,2),'LineStyle', 'none','LineWidth',1,'Color','r','CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
% errorbar([1 2 3 4] + offset, dataStats.mean_area(:,3), dataStats.sem_area(:,3),'LineStyle', 'none','LineWidth',1,'Color',[0.4 0.4 0.4],'CapSize',0,'Marker','.','MarkerSize',meanMarkSize);
% ylabel("Crenation area");
% xticks(1:4);
% xticklabels({'E16-17','P0-2','P6-8','P10-12'});

%% stats for frequency
%E16 and P0 are ns because values are zeros
[p,avova00,stats] = anova1(dataBlock.freqs(dataBlock.Var1==1,:));
c00 = multcompare(stats)

[p,anova0,stats] = anova1(dataBlock.freqs(dataBlock.Var1==1.7,:));
c0 =multcompare(stats)

[p,anova0,stats] = anova1(dataBlock.freqs(dataBlock.Var1==2.3,:));
c =multcompare(stats)


%P7
[p,anova111,stats] = anova1(dataBlock.freqs(dataBlock.Var1==3,:));
c1 =multcompare(stats);
[p,anova222,stats] = anova1(dataBlock.freqs(dataBlock.Var1==4,:));
c2 = multcompare(stats)

%%
[mean, sem] = grpstats(dataBlock.freqs(dataBlock.Var1==2.3,1),ones(6,1),{'mean','sem'})
[mean, sem] = grpstats(dataBlock.freqs(dataBlock.Var1==3,1),ones(size(dataBlock.freqs(dataBlock.Var1==3,1),1),1),{'mean','sem'})
[mean, sem] = grpstats(dataBlock.freqs(dataBlock.Var1==4,1),ones(size(dataBlock.freqs(dataBlock.Var1==4,1),1),1),{'mean','sem'})

%% draw Ca2+ events
%E16
num = 1;
fileList = loadFileList('.\Data\E16\*.tif');
dataFileList = loadFileList('.\Data\E16\*.mat');

img = loadFirstTif(fileList{num},8);
data = load(dataFileList{num});
crens = data.Crens;
[h,h1] = drawCrenations(crens.crenationImg(:,:,crens.locs < 300),crens.locs(crens.locs < 300),300,uint8(img));
figure(h);
print('.\EPS Panels\E16_baseline_crenations.pdf','-dpdf');
figure(h1);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\E16_baseline_crenations_timeline.eps');


[h2,h3] = drawCrenations(crens.crenationImg(:,:,crens.locs > 480),crens.locs(crens.locs > 780)-480,300,uint8(img));
figure(h2);
print('.\EPS Panels\E16_MRS_crenations.pdf','-dpdf');
figure(h3);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\E16_MRS_crenations_timeline.eps');

%%
%P0 Apex
num = 5;
fileList = loadFileList('.\Data\P0-2\*.tif');
dataFileList = loadFileList('.\Data\P0-2\*.mat');

img = loadFirstTif(fileList{num},8);
data = load(dataFileList{num});
crens = data.Crens;
[h,h1] = drawCrenations(crens.crenationImg(:,:,crens.locs < 300),crens.locs(crens.locs < 300),300,uint8(img));
figure(h);
print('.\EPS Panels\P1_baseline_crenations.pdf','-dpdf');
figure(h1);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P1_baseline_crenations_timeline.eps');


[h2,h3] = drawCrenations(crens.crenationImg(:,:,crens.locs > 600 & crens.locs <900),crens.locs(crens.locs > 600 & crens.locs <900)-600,300,uint8(img));
figure(h2);
print('.\EPS Panels\P1_MRS_crenations.pdf','-dpdf');
figure(h3);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P1_MRS_crenations_timeline.eps');
%%
%P0 Base
num = 1;
fileList = loadFileList('.\Data\P0-2\base\*.tif');
dataFileList = loadFileList('.\Data\P0-2\base\*.mat');

img = loadFirstTif(fileList{num},8);
data = load(dataFileList{num});
crens = data.Crens;
[h,h1] = drawCrenations(crens.crenationImg(:,:,crens.locs < 300),crens.locs(crens.locs < 300),300,uint8(img));
figure(h);
print('.\EPS Panels\P1_base_baseline_crenations.pdf','-dpdf');
figure(h1);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P1_base_baseline_crenations_timeline.eps');


[h2,h3] = drawCrenations(crens.crenationImg(:,:,crens.locs > 600 & crens.locs <900),crens.locs(crens.locs > 600 & crens.locs <900)-600,300,uint8(img));
figure(h2);
print('.\EPS Panels\P1_base_MRS_crenations.pdf','-dpdf');
figure(h3);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P1_base_MRS_crenations_timeline.eps');
%%
%P7
num = 2;
fileList = loadFileList('.\Data\P7-8\*.tif');
dataFileList = loadFileList('.\Data\P7-8\*.mat');

img = loadFirstTif(fileList{num},8);
data = load(dataFileList{num});
crens = data.Crens;
[h,h1] = drawCrenations(crens.crenationImg(:,:,crens.locs < 300),crens.locs(crens.locs < 300),300,uint8(img));
figure(h);
print('.\EPS Panels\P7_baseline_crenations.pdf','-dpdf');
figure(h1);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P7_baseline_crenations_timeline.eps');


[h2,h3] = drawCrenations(crens.crenationImg(:,:,crens.locs > 600 & crens.locs <900),crens.locs(crens.locs > 600 & crens.locs <900)-600,300,uint8(img));
figure(h2);
print('.\EPS Panels\P7_MRS_crenations.pdf','-dpdf');
figure(h3);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P7_MRS_crenations_timeline.eps');

%%
%P11
num = 2;
fileList = loadFileList('.\Data\P10-12\*.tif');
dataFileList = loadFileList('.\Data\P10-12\*.mat');

img = loadFirstTif(fileList{num},8);
data = load(dataFileList{num});
crens = data.Crens;
[h,h1] = drawCrenations(crens.crenationImg(:,:,crens.locs < 300),crens.locs(crens.locs < 300),300,uint8(img));
figure(h);
print('.\EPS Panels\P11_baseline_crenations.pdf','-dpdf');
figure(h1);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P11_baseline_crenations_timeline.eps');

img2 = loadFirstTif(fileList{num},8,900);
[h2,h3] = drawCrenations(crens.crenationImg(:,:,crens.locs > 900 & crens.locs <1200),crens.locs(crens.locs > 900 & crens.locs <1200)-900,300,uint8(img));
figure(h2);
print('.\EPS Panels\P11_MRS_crenations.pdf','-dpdf');
figure(h3);
figQuality(gcf,gca,[3 0.5]);
export_fig('.\EPS Panels\P11_MRS_crenations_timeline.eps');



%%

function [freqs, areas, areasmedian] = getStats(crens, conditionTimes)
    freqs = []; areas = []; areasmedian = [];
    for i = 1:size(crens,2)
        for j = 1:size(conditionTimes,1)
            timeMin = (conditionTimes(j,2) - conditionTimes(j,1))/60;
            locs = crens(i).locs;
            locs = locs > conditionTimes(j,1) & locs < conditionTimes(j,2);
            if locs == 0
                locs = [];
            end
            if ~isempty(locs)
                freqs(i,j) = size(crens(i).locs(locs),1)/timeMin ;
                areas(i,j) = nanmean(crens(i).areas(locs),1);
                areasmedian(i,j) = nanmedian(crens(i).areas(locs),1);
            else
                freqs(i,j) = 0;
                areas(i,j) = NaN;
                areasmedian(i,j) = NaN;
            end
        end
    end
end
