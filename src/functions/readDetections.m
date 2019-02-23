function frames = readDetections(detectionDir)

    detectionFile = [detectionDir,'detections.mat'];
    allframedata = load(detectionFile);
    allframedata = allframedata.detections;


    frames = struct([]);
    ifc = iscell(allframedata);
    if ifc 
      numframes = length(allframedata);
    else
      [numframes,dim1,dim2]  = size(allframedata);
      
    end
      
    %disp(size(allframedata));

    for f = 1 : numframes
        
      if ifc
        data = allframedata{f};
      else
        data = reshape(allframedata(f,:,:),dim1,dim2);
      end
    %  scores = [data(:,9:end),data(:,8)];
    %  labels = data(:,2);
      bboxes = data(:,4:7);
      boxes =round([bboxes(:,1)*320  bboxes(:,2)*240 bboxes(:,3)*320 bboxes(:,4)*240]);
    %  indexes = data(:,1);
      frames(f).boxes = boxes;
      frames(f).labels = data(:,2);
      frames(f).indexes = data(:,1);
      frames(f).scores = [data(:,9:end),data(:,8)];
    end

end