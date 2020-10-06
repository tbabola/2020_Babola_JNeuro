defaultDir = 'M:\Projects and Analysis\Papers\Development'; %change if in a different location
cd(defaultDir);
addpath(genpath('MATLAB Functions'));
cd(['MRS2500 over development ISCs']);

%%
%This script analyzes data recorded from ISCs under baseline conditions
%(first 300 s of recording) at all ages

%parameters for analysis
pulse = -100; %test pulse in pA for resistance measurement
pulseTime = [0 0.01 0.2 0.26]; %periods to measure resistance [baseline_start baseline_end ss_start ss_end]
peakProminenceThr = 25; %threshold in pA
conditions = {'Baseline'}; %{'Baseline','MRS','Washout'};
conditionTimes = [0 300]; %[0 300; 540 840; 1800 2100];

%%dir 
dname = '.\Data\WT\*\*\*';
folderlist = dir(dname);
folderlist = folderlist(~contains({folderlist.name},'.'))
analyze = [1 1 1 1 1];

if exist('cell') | exist('cellS') 
   clear 'cell' 'cellS'; 
end

for i = 1:size(folderlist,1)
    if analyze(1) %Input resistance measurement (baseline)
        fileList = loadFileList([folderlist(i).folder '\' folderlist(i).name '\*_rin.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).bl_R = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).bl_R = NaN;
        end
        disp('Rin complete');
    end
        
    if analyze(2) %Input resistance measurement right before experiment
        fileList = loadFileList([folderlist(i).folder '\' folderlist(i).name '\*_p.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).preR = calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).preR = NaN;
        end
    end
    
    if analyze(3) %analysis of gap-free current clamp recording for Vm
        fileList = loadFileList([folderlist(i).folder '\' folderlist(i).name '\*_gfi.abf']);
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
    
    if analyze(4) %analyses of voltage clamp recording (Vhold = -90mV) with peak detection
        fileList = loadFileList([folderlist(i).folder '\' folderlist(i).name '\*_gfv.abf']);
 
        if ~isempty(fileList)
            saveList = fileList;
            saveFile = 1;
            [d,time,SR] = loadPclampData(fileList{1});
            out = msbackadj(time, -d, 'WindowSize', 60,'StepSize',30);
            out = lowpass(out,5,5000);
            %figure for assessing basleine
%              figure; plot(time,out);
%              hold on;
%              plot(time,d-out);
            [pks,locs] = findpeaks(out,time,'MinPeakProminence', peakProminenceThr);
            %uncomment for peaks graph
            %findpeaks(out,time,'MinPeakProminence', peakProminenceThr)
            cellS(i).rawLocPk = [locs pks];
            
            if exist('stats')
                clear 'stats'; 
            end
            for j = 1:size(conditionTimes,1)
                if conditionTimes(j,2) <= time(end)+1
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
            saveFile = 0; 
        end

    end
   
    if analyze(5) %Input resistance measurement obtained after the experiment
        fileList = loadFileList([folderlist(i).folder '\' folderlist(i).name '\*_pp.abf']);
        if ~isempty(fileList)
            [d,time] = loadPclampData(fileList{1});
            cellS(i).postR= calcRm(d,time,pulse,pulseTime);
        else
            cellS(i).postR= NaN;
        end
    end
    
    if saveFile
        tempCell = cellS(i);
        [d,f,ext] = fileparts(saveList{1});
        save([d '\' f '_stats.mat'],'tempCell');
    end
    
    
end
