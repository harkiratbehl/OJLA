function frames = readDetections_f(videoDetDirBase,videoDetDirTop)

    detectionFile = [videoDetDirBase,'detections.mat'];
    baseframedata = load(detectionFile);
    baseframedata = baseframedata.detections;

    detectionFile = [videoDetDirTop,'detections.mat'];
    topframedata = load(detectionFile);
    topframedata = topframedata.detections;

    frames = struct([]);
    ifcB = iscell(baseframedata);
    if ifcB 
      numframesB = length(baseframedata);
    else
      [numframesB,dimB1,dimB2]  = size(baseframedata);
    end

    ifcT = iscell(topframedata);
    if ifcT 
      numframesT = length(topframedata);
    else
      [numframesT,dimT1,dimT2]  = size(topframedata);
    end

    %disp(size(baseframedata))
    %disp(size(topframedata))
    try
        assert(numframesT == numframesB);
    catch
        ru;
    end
    for f = 1 : numframesT
      if ifcB
        dataB = baseframedata{f};
      else
        dataB = reshape(baseframedata(f,:,:),dimB1,dimB2);
      end
      
      if ifcT
        dataT = topframedata{f};
      else
        dataT = reshape(topframedata(f,:,:),dimT1,dimT2);
      end
      
    %  scores = [data(:,9:end),data(:,8)];
    %  labels = data(:,2);
      bboxes = dataB(:,4:7);
      boxesB =round([bboxes(:,1)*320  bboxes(:,2)*240 bboxes(:,3)*320 bboxes(:,4)*240]);
      bboxes = dataT(:,4:7);
      boxesT =round([bboxes(:,1)*320  bboxes(:,2)*240 bboxes(:,3)*320 bboxes(:,4)*240]);
    %  indexes = data(:,1);
      frames(f).baseBoxes = boxesB;
      frames(f).topBoxes = boxesT;
      frames(f).baseScores = [dataB(:,9:end),dataB(:,8)];
      frames(f).topScores = [dataT(:,9:end),dataT(:,8)];

      frames(f).catboxes = [frames(f).baseBoxes;frames(f).topBoxes];
      frames(f).catscores = [frames(f).baseScores;frames(f).topScores];
    end

end