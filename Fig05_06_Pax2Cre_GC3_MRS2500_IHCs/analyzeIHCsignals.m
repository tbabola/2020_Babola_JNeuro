%%load IHC data
list = loadFileList('./Data/P0 MRS2500/*IHC*');

for i = 1:size(list,1)
    l = load(list{i});
    IHCstruct = l.IHCstruct;
    rois = IHCstruct.smoothroi;
    labeled = bwlabel(IHCstruct.thrroi,8);
    for j = 1:max(labeled,[],'all')
        [i1,i2] = ind2sub(size(labeled),find(labeled(:)==j)); %get x,t coordinates for label i
        if max(i2)-min(i2) < 4
            labeled(labeled == j) = 0;
        end
    end
    labeled = bwlabel(labeled,8);
    
    if max(labeled,[],'all') > 0
        for j = 1:max(labeled,[],'all')
            [i1,i2] = ind2sub(size(labeled),find(labeled(:)==j)); %get x,t coordinates for label i
            IHCstruct.event(j).timeStart = min(i1);
            IHCstruct.event(j).timeEnd = max(i1);
            IHCstruct.event(j).eventDuration = max(i1) - min(i1);
            IHCstruct.event(j).area = max(i2) - min(i2);
            IHCstruct.event(j).maxAmplitude = max(rois(find(labeled == j)),[],'all');
        end
    else
        IHCstruct.event = [];
    end
    
    IHCstruct.corr = corr(rois);
    save(list{i}, 'IHCstruct');
    %figure;
    %imagesc([IHCstruct.thrroi labeled>1]);
end
