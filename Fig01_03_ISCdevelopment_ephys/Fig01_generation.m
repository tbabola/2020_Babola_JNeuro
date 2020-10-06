%% Figure 1: ISC Spontaneous Currents
defaultDir = 'D:\Projects and Analysis\Papers\Development'; %change if in a different location
%defaultDir = 'C:\Users\Travis\Desktop\Development Paper';
cd(defaultDir);
addpath(genpath('MATLAB Functions'));
cd(['Fig01_03_ISCdevelopment_ephys']);

%% Panel B
apexFiles = loadFileList('./Data/WT/E16/Apex/*/*_gfv.abf');
baseFiles = loadFileList('./Data/WT/E16/Base/*/*_gfv.abf');
E16_handles = {};
for i = 1:size(apexFiles,1)
    [d, time] = loadPclampData(apexFiles{i});
    dsm = smooth(d,1000);
    dsm = msbackadj(time, -dsm,'WindowSize',20,'StepSize',20);
    figure; E16handles{i} = gcf;
    plot(time,-dsm,'Color','k'); hold on;
    [d, time] = loadPclampData(baseFiles{i});
    dsm = smooth(d,1000);
    dsm = msbackadj(time, -dsm,'WindowSize',20,'StepSize',20);
    plot(time,-dsm - 250,'Color','k');
    ylim([-2500 250]);
    xlim([0 300]);
end

handleTheSubplot(E16handles,[3 1]);
figQuality(gcf,gca,[4,3]);
export_fig('./EPS Panels/E16_Base_Apex_Recordings.eps');
%% Panel C
P0Files = loadFileList('./Data/WT_MRS/P1_2/*/*_gfv.abf');
P7Files = loadFileList('./Data/WT_MRS/P7_8/*/*_gfv.abf');
P11Files = loadFileList('./Data/WT_MRS/P11_12/*/*_gfv.abf');
files = {P0Files{3}; P7Files{1}; P11Files{3}};

handles = {};
for i = 1:size(files,1)
    [d, time] = loadPclampData(files{i});
    dsm = smooth(d,1000);
    dsm = msbackadj(time, -dsm,'WindowSize',20,'StepSize',20);
    figure; handles{i} = gcf;
    plot(time,-dsm,'Color','k'); hold on;
    ylim([-2500 250]);
    xlim([0 300]);
end

handleTheSubplot(handles,[3 1]);
figQuality(gcf,gca,[4,3]);
export_fig('./EPS Panels/Older_Apex_Recordings.eps');
%% Panel D, integral
[apex,base,all] = getAgesLocation(); %load data for all of Panel D

labels = {'E14-16','P0-2','P7-8','P10-12'};
l_green = [192 192 255]/255;
d_green = [0 0 255]/255;
l_mag = [208 150 150]/255;
d_mag = [145 22 26]/255;

h = figure; scatter(apex.Int.Group,apex.Int.Int,5,l_green,'filled','jitter','on','jitterAmount',0.15); hold on;
scatter(base.Int.Group,base.Int.Int,5,l_mag,'filled','jitter','on','jitterAmount',0.15); 
intStats = grpstats(apex.Int,'Group',{'mean','sem'}); 
errorbar(intStats.Group, intStats.mean_Int, intStats.sem_Int,'LineStyle', 'none','LineWidth',1,'Color',d_green,'CapSize',0,'Marker','.','MarkerSize',12);
baseintStats = grpstats(base.Int,'Group',{'mean','sem'}); 
errorbar(baseintStats.Group, baseintStats.mean_Int, baseintStats.sem_Int,'LineStyle', 'none','LineWidth',1,'Color',d_mag,'CapSize',0,'Marker','.','MarkerSize',12);

xticks(1:4);
xlim([0 4.5]);
xticklabels(labels);
xtickangle(45);
figQuality(gcf,gca,[3 2]);
ylabel('Charge transfer (-pC)');
%% Panel D, frequency
h1 = figure; scatter(apex.Freq.Group,apex.Freq.Freq,5,l_green,'filled','jitter','on','jitterAmount',0.15); hold on;
scatter(base.Freq.Group,base.Freq.Freq,5,l_mag,'filled','jitter','on','jitterAmount',0.15); 
freqStats = grpstats(apex.Freq,'Group',{'mean','sem'}); 
errorbar(freqStats.Group, freqStats.mean_Freq, freqStats.sem_Freq,'LineStyle', 'none','LineWidth',1,'Color',d_green,'CapSize',0,'Marker','.','MarkerSize',12);
basefreqStats = grpstats(base.Freq,'Group',{'mean','sem'}); 
errorbar(basefreqStats.Group, basefreqStats.mean_Freq, basefreqStats.sem_Freq,'LineStyle', 'none','LineWidth',1,'Color',d_mag,'CapSize',0,'Marker','.','MarkerSize',12);

xticks(1:4);
xlim([0 4.5]);
xticklabels(labels);
xtickangle(45);
figQuality(gcf,gca,[3 2]);
ylabel('Events per minute');
%% Panel D, amplitude
h2 = figure; scatter(apex.Amp.Group,apex.Amp.Amp,5,l_green,'filled','jitter','on','jitterAmount',0.15); hold on;
scatter(base.Amp.Group,base.Amp.Amp,5,l_mag,'filled','jitter','on','jitterAmount',0.15); 
ampApexStats = grpstats(apex.Amp,'Group',{'mean','sem'});
errorbar(ampApexStats.Group, ampApexStats.mean_Amp, ampApexStats.sem_Amp,'LineStyle', 'none','LineWidth',1,'Color',d_green,'CapSize',0,'Marker','.','MarkerSize',12);

ampBaseStats = grpstats(base.Amp,'Group',{'mean','sem'});
errorbar(1.2, ampBaseStats.mean_Amp(1), ampBaseStats.sem_Amp(1),'LineStyle', 'none','LineWidth',1,'Color',d_mag,'CapSize',0,'Marker','.','MarkerSize',12);

xticks(1:4);
xlim([0 4.5]);
xticklabels(labels);
xtickangle(45);
ylabel('Amplitude (-pA)');
figQuality(gcf,gca,[3 2]);
%% Panel D, input resistance
h4 = figure; scatter(apex.Rin.Group,apex.Rin.Rin,5,l_green,'filled','jitter','on','jitterAmount',0.15); hold on;
scatter(base.Rin.Group,base.Rin.Rin,5,l_mag,'filled','jitter','on','jitterAmount',0.15); 
ampStats = grpstats(apex.Rin,'Group',{'mean','sem'}); hold on;
errorbar(ampStats.Group, ampStats.mean_Rin, ampStats.sem_Rin,'LineStyle', 'none','LineWidth',1,'Color',d_green,'CapSize',0,'Marker','.','MarkerSize',12);
ampBaseStats = grpstats(base.Rin,'Group',{'mean','sem'});
errorbar(1.2, ampBaseStats.mean_Rin(1), ampBaseStats.sem_Rin(1),'LineStyle', 'none','LineWidth',1,'Color',d_mag,'CapSize',0,'Marker','.','MarkerSize',12);

xticks(1:4);
xlim([0 4.5]);
xticklabels(labels);
xtickangle(45);
ylabel('Input resistance (M\Omega)');
figQuality(gcf,gca,[3 2]);
 handleTheSubplot({h1,h2,h,h4},[2 2])
 figQuality(gcf,gca,[3.75 3]);
 export_fig('.\EPS Panels\Spont_Act_Onset_with_Base.eps');

%% statsistics
%comparison of base to apex at E16, data are clearly non-normal so ranksum was used
[h,p,stats] = ranksum(apex.Freq.Freq(apex.Freq.Group == 0.8),base.Freq.Freq(base.Freq.Group == 1.2))
[h,p,stats] = ranksum(apex.Amp.Amp(apex.Amp.Group == 0.8),base.Amp.Amp(base.Amp.Group == 1.2))
[h,p,stats] = ranksum(apex.Int.Int(apex.Int.Group == 0.8),base.Int.Int(base.Int.Group == 1.2))
[h,p,stats] = ranksum(apex.Rin.Rin(apex.Rin.Group == 0.8),base.Rin.Rin(base.Rin.Group == 1.2))

%comparison of activity across the apical recordings
tempGroup = apex.Freq.Group;
[p,~,stats] = anova1(apex.Freq.Freq,tempGroup);
c = multcompare(stats);

[p,~,stats] = anova1(apex.Amp.Amp,tempGroup);
c1 = multcompare(stats);

[p,~,stats] = anova1(apex.Int.Int,tempGroup);
c2 = multcompare(stats);

[p,~,stats] = anova1(apex.Rin.Rin,tempGroup);
c3 = multcompare(stats);

%% Print stats
display("Frequencies")
[statarray, sem] = grpstats(apex.Freq.Freq,tempGroup,{'mean','sem'})
display("Amplitudes")
[statarray, sem] = grpstats(apex.Amp.Amp,tempGroup,{'mean','sem'})
display("Integrals")
[statarray, sem] = grpstats(apex.Int.Int,tempGroup,{'mean','sem'})
%%

function [apex,base,all] = getAgesLocation()
    %apical
    cells_E14 = loadCellStructs('.\Data\WT\E14\Apex\*\*_gfv_stats.mat');
    cells_E15 = loadCellStructs('.\Data\WT\E15\*\*\*_gfv_stats.mat');
    cells_E16 = loadCellStructs('.\Data\WT\E16\Apex\*\*_gfv_stats.mat');
    %cells_E16 = [cells_E16 loadCellStructs('.\Data\WT_MRS\E16_17\200123*\*_gfv_stats.mat')];
    %cells_E17 = loadCellStructs('.\Data\WT_MRS\E16_17\200116*\*_gfv_stats.mat');
    cells_P1 = loadCellStructs('.\Data\WT_MRS\P1_2\*\*_gfv_stats.mat');
    cells_P7 = loadCellStructs('.\Data\WT_MRS\P7_8\*\*_gfv_stats.mat');
    cells_P11 = loadCellStructs('.\Data\WT_MRS\P11_12\*\*_gfv_stats.mat');

    bigInt = [getIntegralBlock(cells_E14,0.8); getIntegralBlock(cells_E16,0.8); ...
                getIntegralBlock(cells_P1,2); getIntegralBlock(cells_P7,3); ...
                getIntegralBlock(cells_P11,4)];
    apexInt = table(bigInt(:,1),bigInt(:,2),'VariableNames',{'Group','Int'});

    bigFreq = [getFreqBlock(cells_E14,0.8); getFreqBlock(cells_E16,0.8); ...
                getFreqBlock(cells_P1,2); getFreqBlock(cells_P7,3); ...
                getFreqBlock(cells_P11,4)];
    apexFreq = table(bigFreq(:,1),bigFreq(:,2),'VariableNames',{'Group','Freq'});

    bigAmp = [getAmpBlock(cells_E14,0.8); getAmpBlock(cells_E16,0.8); ...
                getAmpBlock(cells_P1,2); getAmpBlock(cells_P7,3); ...
                getAmpBlock(cells_P11,4)];
    apexAmp = table(bigAmp(:,1), bigAmp(:,2), 'VariableNames',{'Group','Amp'});
    
    bigRin = [getRinBlock(cells_E14,0.8); getRinBlock(cells_E16,0.8); ...
                getRinBlock(cells_P1,2); getRinBlock(cells_P7,3); ...
                getRinBlock(cells_P11,4)];
    apexRin = table(bigRin(:,1), bigRin(:,2), 'VariableNames',{'Group','Rin'});
    
    
    apex = struct();
    apex.Int = apexInt;
    apex.Freq = apexFreq;
    apex.Amp = apexAmp;
    apex.Rin = apexRin;
    
    %base
    cells_E14 = loadCellStructs('.\Data\WT\E14\Base\*\*_gfv_stats.mat');
    cells_E16 = loadCellStructs('.\Data\WT\E16\Base\*\*_gfv_stats.mat');
    bigInt = [getIntegralBlock(cells_E14,1.2);  getIntegralBlock(cells_E16,1.2)];
    baseInt = table(bigInt(:,1),bigInt(:,2),'VariableNames',{'Group','Int'});

    bigFreq = [getFreqBlock(cells_E14,1.2); getFreqBlock(cells_E16,1.2)];
    baseFreq = table(bigFreq(:,1),bigFreq(:,2),'VariableNames',{'Group','Freq'});

    bigAmp = [getAmpBlock(cells_E14,1.2); getAmpBlock(cells_E16,1.2)];
    baseAmp = table(bigAmp(:,1), bigAmp(:,2), 'VariableNames',{'Group','Amp'});
    
    bigRin = [getRinBlock(cells_E14,1.2); getRinBlock(cells_E16,1.2)];
    baseRin = table(bigRin(:,1), bigRin(:,2), 'VariableNames',{'Group','Rin'});
    
    base = struct();
    base.Int = baseInt;
    base.Freq = baseFreq;
    base.Amp = baseAmp;
    base.Rin = baseRin;
    
    
    %all
    cells_E14 = loadCellStructs('.\Data\WT\E14\*\*\*_gfv_stats.mat');
    cells_E16 = loadCellStructs('.\Data\WT\E16\*\*\*_gfv_stats.mat');

    bigInt = [getIntegralBlock(cells_E14,1);  getIntegralBlock(cells_E16,1); ...
                getIntegralBlock(cells_P1,2); getIntegralBlock(cells_P7,3); ...
                getIntegralBlock(cells_P11,4)];
    allInt = table(bigInt(:,1),bigInt(:,2),'VariableNames',{'Group','Int'});

    bigFreq = [getFreqBlock(cells_E14,1);  getFreqBlock(cells_P7,3); ...
                getFreqBlock(cells_P11,4)];
    allFreq = table(bigFreq(:,1),bigFreq(:,2),'VariableNames',{'Group','Freq'});

    bigAmp = [getAmpBlock(cells_E14,1);  getAmpBlock(cells_E16,1); ...
                getAmpBlock(cells_P1,2); getAmpBlock(cells_P7,3); ...
                getAmpBlock(cells_P11,4)];
    allAmp = table(bigAmp(:,1), bigAmp(:,2), 'VariableNames',{'Group','Amp'});
    
    bigRin = [getRinBlock(cells_E14,1);  getRinBlock(cells_E16,1); ...
                getRinBlock(cells_P1,2); getRinBlock(cells_P7,3); ...
                getRinBlock(cells_P11,4)];
    allRin = table(bigRin(:,1), bigRin(:,2), 'VariableNames',{'Group','Rin'});
    
    
    all = struct();
    all.Int = allInt;
    all.Freq = allFreq;
    all.Amp = allAmp;
    all.Rin = allRin;
    
    
end

function integrals = getIntegralBlock(cells,ageIndex)
    for i = 1:size(cells,2)
        integrals(i,:) = [ageIndex, cells(i).stats(1).integral/300];
    end
end

function frequencies = getFreqBlock(cells,ageIndex)
    for i = 1:size(cells,2)
        frequencies(i,:) = [ageIndex, cells(i).stats(1).numPeaks/5];
    end
end

function amps = getAmpBlock(cells,ageIndex)
    for i = 1:size(cells,2)
        amps(i,:) = [ageIndex, cells(i).stats(1).meanAmp];
    end
end

function rins = getRinBlock(cells,ageIndex)
    for i = 1:size(cells,2)
        rins(i,:) = [ageIndex, cells(i).preR];
    end
end