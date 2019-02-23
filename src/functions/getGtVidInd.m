function gtVidInd = getGtVidInd(video,videoName)
    gtVidInd = -1;
    for i=1:length(video)
        vidid = video(i).name;
        vidid = strsplit(vidid,'/');
        vidid = vidid{2};
        if strcmp(vidid,videoName)
            gtVidInd = i;
            break;
        end
    end
end