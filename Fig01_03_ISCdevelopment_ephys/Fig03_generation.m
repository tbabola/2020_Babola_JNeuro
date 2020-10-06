%% Figure 3: MRS2500
defaultDir = 'D:\Projects and Analysis\Papers\Development'; %change if in a different location
%defaultDir = 'C:\Users\Travis\Desktop\Development Paper'
cd(defaultDir);
addpath(genpath('MATLAB Functions'));
cd(['Fig01_03_ISCdevelopment_ephys']);

%% Panel A: Example traces
fileList = loadFileList('.\Data\WT_MRS\E16_17\*\*gfv.abf');

[d,time0,Fs] = loadPclampData(fileList{5}); %11
dsm = smooth(d,100);
dsm16_5 = msbackadj(time0, -dsm,'WindowSize',20);


fileList = loadFileList('.\Data\WT_MRS\P1_2\*\*gfv.abf');

[d,time1,Fs] = loadPclampData(fileList{11}); %11
dsm = smooth(d,100);
dsm1_2 = msbackadj(time1, -dsm,'WindowSize',20);


fileList = loadFileList('.\Data\WT_MRS\P7_8\*\*gfv.abf');
[d,time2,Fs] = loadPclampData(fileList{2});
dsm = smooth(d,100);
dsm7_8 = msbackadj(time2, -dsm,'WindowSize',20);

fileList = loadFileList('.\Data\WT_MRS\P11_12\*\*gfv.abf');
[d,time3,Fs] = loadPclampData(fileList{10});
dsm = smooth(d,100);
dsm11_12 = msbackadj(time3, -dsm,'WindowSize',20);
%% plot the figure
spacer = 1000;
figure(1); hold off; plot(time0+100,-dsm16_5, 'Color','k'); hold on; 
plot(time1,-dsm1_2 - spacer, 'Color','k');
plot(time2,-dsm7_8 - 2*spacer,'Color','k');
plot(time3,-dsm11_12 - 3 * spacer,'Color','k');
xlim([200 780]);
figQuality(gcf,gca,[3 4]);
export_fig('.\EPS Panels\A_exampleTraces.eps')

%% Panel B Frequency, amplitude, and integral plots
% = loadFileList(['.\Data\WT_MRS\P1_2\*\*_gfv_stats.mat']);
cellS16_5 = loadCellStructs('.\Data\WT_MRS\E16_17\*\*_gfv_stats.mat');
cellS1_2 = loadCellStructs('.\Data\WT_MRS\P1_2\*\*_gfv_stats.mat');
cellS7_8 = loadCellStructs('.\Data\WT_MRS\P7_8\*\*_gfv_stats.mat');
cellS11_12 = loadCellStructs('.\Data\WT_MRS\P11_12\*\*_gfv_stats.mat')

[pkstatBlock16_5, ampstatBlock16_5,intstatBlock16_5, rinstatblock16_5] = getstatBlock(cellS16_5);
[pkstatBlock1_2, ampstatBlock1_2,intstatBlock1_2, rinstatblock1_2] = getstatBlock(cellS1_2);
[pkstatBlock7_8, ampstatBlock7_8,intstatBlock7_8,rinstatblock7_8] = getstatBlock(cellS7_8);
[pkstatBlock11_12, ampstatBlock11_12,intstatBlock11_12, rinstatblock11_12] = getstatBlock(cellS11_12);

pkstatBlock(6:7,3) = NaN;
ampstatBlock(6:7,3) = NaN;
intstatBlock(6:7,3) = NaN;

h1 = compareDev(pkstatBlock16_5,pkstatBlock1_2,pkstatBlock7_8, pkstatBlock11_12,{'E15-16','P0-2','P6-8','P10-12'}, 'Events per minute',[1.5 0.9],[12 10]);
figure(h1); 
ylim([0 50]);
yticks(0:25:50);
%export_fig('.\EPS Panels\B_Frequency.eps');

h2 = compareDev(ampstatBlock16_5,ampstatBlock1_2,ampstatBlock7_8, ampstatBlock11_12,{'E15-16','P0-2','P6-8','P10-12'}, 'Amplitude (-pA)',[1.5 0.9],[12 10]);
figure(h2); 
ylim([0 1500]);
yticks(0:500:1500);
%export_fig('.\EPS Panels\C_Amplitude.eps');

h3 = compareDev(intstatBlock16_5,intstatBlock1_2,intstatBlock7_8, intstatBlock11_12,{'E15-16','P0-2','P6-8','P10-12'}, 'Integral (-pC)',[1.5 0.9],[12 10]);
figure(h3); 
ylim([0 1000]);
yticks(0:500:1000);
%export_fig('.\EPS Panels\C_Amplitude.eps');

h4 = handleTheSubplot({h1,h2,h3},[3 1]);
figQuality(gcf,gca,[2 3]);
export_fig('.\EPS Panels\B_D_Freq_Amp_Integral.eps');


%% statistics
%t-test with Bonferroni, we don't care about comparisons between groups
%really

[~,p,~,stats] = ttest(pkstatBlock16_5(:,1),pkstatBlock16_5(:,2))

[~,p,~,stats] = ttest(ampstatBlock16_5(:,1),ampstatBlock16_5(:,2))

[~,p,~,stats] = ttest(intstatBlock16_5(:,1),intstatBlock16_5(:,2))


disp('P0')

[~,p,~,stats]  = ttest(pkstatBlock1_2(:,1),pkstatBlock1_2(:,2))
[~,p,~,stats]  = ttest(ampstatBlock1_2(:,1),ampstatBlock1_2(:,2))
[~,p,~,stats] = ttest(intstatBlock1_2(:,1),intstatBlock1_2(:,2))

disp('P7')
[~,p,~,stats]  = ttest(pkstatBlock7_8(:,1),pkstatBlock7_8(:,2))
[~,p,~,stats]  = ttest(ampstatBlock7_8(:,1),ampstatBlock7_8(:,2))
[~,p,~,stats]  = ttest(intstatBlock7_8(:,1),intstatBlock7_8(:,2))

disp('P11')
[~,p,~,stats]  = ttest(pkstatBlock11_12(:,1),pkstatBlock11_12(:,2))
[~,p,~,stats]  = ttest(ampstatBlock11_12(:,1),ampstatBlock11_12(:,2))
[~,p,~,stats]  = ttest(intstatBlock11_12(:,1),intstatBlock11_12(:,2))
%%
%%print stats
display("E16.5");
display("Frequencies");
tempGroup = [ones(size(pkstatBlock16_5,1),1); 2*ones(size(pkstatBlock16_5,1),1)]
[statarray, sem] = grpstats([pkstatBlock16_5(:,1);pkstatBlock16_5(:,2)],tempGroup,{'mean','sem'})
display("Amplitudes")
[statarray, sem] = grpstats([ampstatBlock16_5(:,1);ampstatBlock16_5(:,2)],tempGroup,{'mean','sem'})
display("Integrals")
[statarray, sem] = grpstats([intstatBlock16_5(:,1);intstatBlock16_5(:,2)],tempGroup,{'mean','sem'})

%%
function [pkstatBlock, ampstatBlock, intstatBlock, rinstatblock] = getstatBlock(cellS)
    for i = 1:size(cellS,2)
        pkstatBlock(i,:) = [cellS(i).stats.numPeaks]/5;
        ampstatBlock(i,:) = [cellS(i).stats.meanAmp];
        intstatBlock(i,:) = [cellS(i).stats.integral]/300;
        rinstatblock(i,:) = [cellS(i).bl_R];
    end
end