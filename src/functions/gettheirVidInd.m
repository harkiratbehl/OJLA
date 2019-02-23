function gtVidInd = gettheirVidInd(video,videoName)
    gtVidInd = -1;
    for i=1:length(video)
        vidid = video(i).videoName;
        if strcmp(vidid,videoName)
            gtVidInd = i;
            break;
        end
    end
end