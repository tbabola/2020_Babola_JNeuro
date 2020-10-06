directory = '.\Data\P0 Base\*ISCdata.mat';
addISCmeanAmp(directory);

directory = '.\Data\P0 Apex Middle\*ISCdata.mat';
addISCmeanAmp(directory);

directory = '.\Data\P0 Apex Apical\*ISCdata.mat';
addISCmeanAmp(directory);


%%
function addISCmeanAmp(directory)
    ISCstructs = loadCellStructs(directory)
    filelist = loadFileList(directory)

    for i = 1:size(ISCstructs,2)
        ISCstruct = ISCstructs(i);
        ISCnum = max(ISCstruct.labelRoi,[],'all')
        labelRoi = ISCstruct.labelRoi;
        for j = 1:ISCnum
            [r,c] = find(labelRoi == j);
            r = unique(r);
            c = unique(c);
            meanSignal = mean(ISCstruct.rois(min(r):max(r),c),2);
            %figure; plot(meanSignal);
            ISCstruct.event(j).meanAmplitude = max(meanSignal);
        end
        ISCstructs(i) = ISCstruct;
        [fp,name,ext] = fileparts(filelist{i})
        newFolder = [fp '\newISCstructs\']
        if ~exist(newFolder, 'dir')
           mkdir(newFolder)
        end
        save([newFolder name '_ISCdata.mat'],'ISCstruct');
    end
end



