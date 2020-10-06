defaultDir = 'M:\Projects and Analysis\Papers\Development'; %change if in a different location
cd(defaultDir);
addpath(genpath('MATLAB Functions'));
cd(['MRS2500 over development ISCs']);
%%
%This script analyzes data recorded from ISCs under baseline and MRS conditions
% in the directory specified by dname and at age specified. Slightly
% different conditions across ages is accounted for by adjusting the
% timing.

%parameters for analysis
pulse = -100; %test pulse in pA for resistance measurement
pulseTime = [0 0.01 0.2 0.26]; %periods to measure resistance [baseline_start baseline_end ss_start ss_end]
peakProminenceThr = 25; %threshold in pA
conditions = {'Baseline','MRS','Washout'};
conditionTimes = [0 300; 540 840; 1800 2100];

%%dir 
dname = '.\Data\WT_MRS\P1_2';
folderlist = dir(dname);
folderlist = folderlist(3:end);
analyze = [1 1 1 1 1];

if exist('cell') | exist('cellS') 
   clear 'cell' 'cellS'; 
end

for i = 1:size(folderlist,1)
    if analyze(1)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_rin.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).bl_R = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).bl_R = NaN;
        end
    end
        
    if analyze(2)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_p.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).preR = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).preR = NaN;
        end
    end
    
    if analyze(3)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_gfi.abf']);
        if ~isempty(fileList)
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, d, 'WindowSize', 60,'StepSize',60);
            cellS(i).Vm = mean(d - out);
        else
            cellS(i).Vm = NaN;
        end
        %%figure for assessing basleine
        % figure; plot(time,d);
        %hold on; plot(time,out);
       % plot(time,d-out);
    end
    
    if analyze(4) %change in inward current
        fileList = loadFileList([dname '\' folderlist(i).name '\*_gfv.abf']);
        saveList = fileList;
        if ~isempty(fileList)
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, -d, 'WindowSize', 60,'StepSize',30);
            out = lowpass(out,5,5000);
            %figure for assessing basleine
            % figure; plot(time,out);
    %         hold on;
    %         plot(time,d-out);
            [pks,locs] = findpeaks(out,time,'MinPeakProminence', peakProminenceThr);
            %uncomment for peaks graph
            %findpeaks(out,time,'MinPeakProminence', peakProminenceThr)
            cellS(i).rawLocPk = [locs pks];
            
            if exist('stats')
                clear 'stats'; 
            end
            for j = 1:size(conditionTimes,1)
                if conditionTimes(j,2) < time(end)
                    stats(j).condition = conditions{j};
                    stats(j).numPeaks = size(pks(locs > conditionTimes(j,1) & locs < conditionTimes(j,2)),1);
                    stats(j).meanAmp = mean(pks(locs > conditionTimes(j,1) & locs < conditionTimes(j,2)),1);
                    cumFunc = cumtrapz(time(time > conditionTimes(j,1) & time < conditionTimes(j,2)), ... 
                                             out(time > conditionTimes(j,1) & time < conditionTimes(j,2)));
                    stats(j).integral = cumFunc(end);
                else
                    stats(j).condition = conditions{j};
                    stats(j).numPeaks = NaN;
                    stats(j).meanAmp = NaN;
                    stats(j).integral = NaN;
                end
                
            end
            
            cellS(i).stats = stats;
        else
            cellS(i).rawLocPk = NaN;
            cellS(i).stats = NaN;
        end

    end
   
    if analyze(5)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_pp.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).postR= calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).postR= NaN;
        end
    end
    
    tempCell = cellS(i);
    [d,f,ext] = fileparts(saveList{1});
    save([d '\' f '_stats.mat'],'tempCell');
    
    
end

%%
pkstatBlock = [];
ampstatBlock = [];
intstatBlock = [];

for i = 1:size(cellS,2)
    pkstatBlock(i,:) = [cellS(i).stats.numPeaks]/5;
    ampstatBlock(i,:) = [cellS(i).stats.meanAmp];
    intstatBlock(i,:) = [cellS(i).stats.integral]/300;
end

dim = [1.6 1.6];
compare3P(pkstatBlock(:,1),pkstatBlock(:,2),pkstatBlock(:,3),conditions,'Events Per Minute',dim,[5 12]);
ylim([0 25]);
yticks(0:5:25);
export_fig('.\EPS Panels\MRSfreq.eps');

compare3P(ampstatBlock(:,1),ampstatBlock(:,2),ampstatBlock(:,3),conditions,'Mean Amp (pA)',dim,[5 12]);
ylim([0 600]);
yticks(0:200:600);
export_fig('.\EPS Panels\MRSamp.eps');

compare3P(intstatBlock(:,1),intstatBlock(:,2),intstatBlock(:,3),conditions,'Charge Transfer (-pA*s)',dim,[5 12]);
ylim([0 400]);
yticks(0:100:400);
export_fig('.\EPS Panels\MRSint.eps');
%%
[p,~,stats] = anova1(pkstatBlock);
[c] = multcompare(stats,'display','off')

[p,~,stats] = anova1(ampstatBlock);
[c1] = multcompare(stats,'display','off')

[p,~,stats] = anova1(intstatBlock);
[c2] = multcompare(stats,'display','off')

%%
%parameters for analysis
pulse = -100; %test pulse in pA for resistance measurement
pulseTime = [0 0.01 0.2 0.26]; %periods to measure resistance [baseline_start baseline_end ss_start ss_end]
peakProminenceThr = 25; %threshold in pA
conditions = {'Baseline','MRS','Washout'};
conditionTimes = [0 300; 540 840; 1800 2100];

%%dir 
dname = '.\Data\WT_MRS\P11_12';
folderlist = dir(dname);
folderlist = folderlist(3:end);
analyze = [1 1 1 1 1];

if exist('cell') | exist('cellS') 
   clear 'cell' 'cellS'; 
end

for i = 1:size(folderlist,1)
    if analyze(1)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_rin.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).bl_R = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).bl_R = NaN;
        end
    end
        
    if analyze(2)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_p.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).preR = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).preR = NaN;
        end
    end
    
    if analyze(3)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_gfi.abf']);
        if ~isempty(fileList)
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, d, 'WindowSize', 60,'StepSize',60);
            cellS(i).Vm = mean(d - out);
        else
            cellS(i).Vm = NaN;
        end
        %%figure for assessing basleine
        % figure; plot(time,d);
        %hold on; plot(time,out);
       % plot(time,d-out);
    end
    
    if analyze(4) %change in inward current
        fileList = loadFileList([dname '\' folderlist(i).name '\*_gfv.abf']);
        saveList = fileList;
        if ~isempty(fileList)
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, -d, 'WindowSize', 60,'StepSize',30);
            out = lowpass(out,5,5000);
            %figure for assessing basleine
            % figure; plot(time,out);
    %         hold on;
    %         plot(time,d-out);
            [pks,locs] = findpeaks(out,time,'MinPeakProminence', peakProminenceThr);
            %uncomment for peaks graph
            %findpeaks(out,time,'MinPeakProminence', peakProminenceThr)
            cellS(i).rawLocPk = [locs pks];
            
            if exist('stats')
                clear 'stats'; 
            end
            for j = 1:size(conditionTimes,1)
                if conditionTimes(j,2) < time(end)
                    stats(j).condition = conditions{j};
                    stats(j).numPeaks = size(pks(locs > conditionTimes(j,1) & locs < conditionTimes(j,2)),1);
                    stats(j).meanAmp = mean(pks(locs > conditionTimes(j,1) & locs < conditionTimes(j,2)),1);
                    cumFunc = cumtrapz(time(time > conditionTimes(j,1) & time < conditionTimes(j,2)), ... 
                                             out(time > conditionTimes(j,1) & time < conditionTimes(j,2)));
                    stats(j).integral = cumFunc(end);
                else
                    stats(j).condition = conditions{j};
                    stats(j).numPeaks = NaN;
                    stats(j).meanAmp = NaN;
                    stats(j).integral = NaN;
                end
                
            end
            
            cellS(i).stats = stats;
        else
            cellS(i).rawLocPk = NaN;
            cellS(i).stats = NaN;
        end

    end
   
    if analyze(5)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_pp.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).postR= calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).postR= NaN;
        end
    end
    
    tempCell = cellS(i);
    [d,f,ext] = fileparts(saveList{1});
    save([d '\' f '_stats.mat'],'tempCell');
    
    
end

%%
pkstatBlock = [];
ampstatBlock = [];
intstatBlock = [];

for i = 1:size(cellS,2)
    pkstatBlock(i,:) = [cellS(i).stats.numPeaks]/5;
    ampstatBlock(i,:) = [cellS(i).stats.meanAmp];
    intstatBlock(i,:) = [cellS(i).stats.integral]/300;
end

dim = [1.6 1.6];
compare3P(pkstatBlock(:,1),pkstatBlock(:,2),pkstatBlock(:,3),conditions,'Events Per Minute',dim,[5 12]);
ylim([0 25]);
yticks(0:5:25);
export_fig('.\EPS Panels\10_12_MRSfreq.eps');

compare3P(ampstatBlock(:,1),ampstatBlock(:,2),ampstatBlock(:,3),conditions,'Mean Amp (pA)',dim,[5 12]);
ylim([0 600]);
yticks(0:200:600);
export_fig('.\EPS Panels\10_12_MRSamp.eps');

compare3P(intstatBlock(:,1),intstatBlock(:,2),intstatBlock(:,3),conditions,'Charge Transfer (-pA*s)',dim,[5 12]);
ylim([0 400]);
yticks(0:100:400);
export_fig('.\EPS Panels\10_12_MRSint.eps');
%%
[p,~,stats] = anova1(pkstatBlock(:,1:2));
[c] = multcompare(stats,'display','off')

[p,~,stats] = anova1(ampstatBlock(:,1:2));
[c1] = multcompare(stats,'display','off')

[p,~,stats] = anova1(intstatBlock(:,1:2));
[c2] = multcompare(stats,'display','off')

%%

%parameters for analysis, P7_8 animals
pulse = -100; %test pulse in pA for resistance measurement
pulseTime = [0 0.01 0.2 0.26]; %periods to measure resistance [baseline_start baseline_end ss_start ss_end]
peakProminenceThr = 25; %threshold in pA
conditions = {'Baseline','MRS','MRS+SP','Washout'};
conditionTimes = [0 300; 540 840; 960 1260; 2100 2400];
%%dir 
dname = '.\Data\WT_MRS\P7_8';
folderlist = dir(dname);
folderlist = folderlist(3:end);
analyze = [1 1 1 1 1];

if exist('cell') | exist('cellS') 
   clear 'cell' 'cellS'; 
end

for i = 1:size(folderlist,1)
    if analyze(1)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_rin.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).bl_R = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).bl_R = NaN;
        end
    end
        
    if analyze(2)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_p.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).preR = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).preR = NaN;
        end
    end
    
    if analyze(3)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_gfi.abf']);
        if ~isempty(fileList)
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, d, 'WindowSize', 60,'StepSize',60);
            cellS(i).Vm = mean(d - out);
        else
            cellS(i).Vm = NaN;
        end
        %%figure for assessing basleine
        % figure; plot(time,d);
        %hold on; plot(time,out);
       % plot(time,d-out);
    end
    
    if analyze(4) %change in inward current
        fileList = loadFileList([dname '\' folderlist(i).name '\*_gfv.abf']);
        saveList = fileList;
        if ~isempty(fileList)
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, -d, 'WindowSize', 60,'StepSize',30);
            out = lowpass(out,5,5000);
            %figure for assessing basleine
            % figure; plot(time,out);
    %         hold on;
    %         plot(time,d-out);
            [pks,locs] = findpeaks(out,time,'MinPeakProminence', peakProminenceThr);
            %uncomment for peaks graph
            %findpeaks(out,time,'MinPeakProminence', peakProminenceThr)
            cellS(i).rawLocPk = [locs pks];
            
            if exist('stats')
                clear 'stats'; 
            end
            for j = 1:size(conditionTimes,1)
                if conditionTimes(j,2) < time(end)
                    stats(j).condition = conditions{j};
                    stats(j).numPeaks = size(pks(locs > conditionTimes(j,1) & locs < conditionTimes(j,2)),1);
                    stats(j).meanAmp = mean(pks(locs > conditionTimes(j,1) & locs < conditionTimes(j,2)),1);
                    cumFunc = cumtrapz(time(time > conditionTimes(j,1) & time < conditionTimes(j,2)), ... 
                                             out(time > conditionTimes(j,1) & time < conditionTimes(j,2)));
                    stats(j).integral = cumFunc(end);
                else
                    stats(j).condition = conditions{j};
                    stats(j).numPeaks = NaN;
                    stats(j).meanAmp = NaN;
                    stats(j).integral = NaN;
                end
                
            end
            
            cellS(i).stats = stats;
        else
            cellS(i).rawLocPk = NaN;
            cellS(i).stats = NaN;
        end

    end
   
    if analyze(5)
        fileList = loadFileList([dname '\' folderlist(i).name '\*_pp.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).postR= calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).postR= NaN;
        end
    end
    
    tempCell = cellS(i);
    [d,f,ext] = fileparts(saveList{1});
    save([d '\' f '_stats.mat'],'tempCell');
    
    
end