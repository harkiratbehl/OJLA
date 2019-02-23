function main()
    
    close all;
    clc;
    clear all;
    clear mex;
    clear is_valid_handle; % to clear init_key
    fileIDc = fopen('src/code_base_path.txt','r');
    codebase_path = fscanf(fileIDc,'%s');
    run([codebase_path '/' 'src/startup']);
    addpath(genpath(fullfile(codebase_path, 'src', 'functions')));
    fclose(fileIDc);
    img_path = [codebase_path , '/datasets/UCF101/images'];
    st_vid_list = load([codebase_path , '/datasets/UCF101/ucf101_annot/testlist.mat']); % video
    videos = st_vid_list.testlist;
    clear st_vid_list;
    st_annot = load([codebase_path , '/datasets/UCF101/ucf101_annot/finalAnnots.mat']); % videos - new modifed APT annot
    Videos = st_annot.annot; clear st_annot;
    num_videos = length(videos);
    num = 0;
    num_actions = 24;
    actions = {'Basketball', 'BasketballDunk', 'Biking', 'CliffDiving', 'CricketBowling', 'Diving', 'Fencing',...
    'FloorGymnastics', 'GolfSwing', 'HorseRiding', 'IceDancing', 'LongJump', 'PoleVault', 'RopeClimbing', 'SalsaSpin', ...
    'SkateBoarding', 'Skiing', 'Skijet', 'SoccerJuggling', 'Surfing', 'TennisSwing', 'TrampolineJumping', 'VolleyballSpiking', 'WalkingWithDog'};
    min_num_frames = {0,20, 18, 15, 19,17,17,17, 18,18,18,17, 29, 17, 18, 17, 21,17, 18,18, 15, 17, 18, 17, 21};
    nms_th = 0.1;
    topk_boxes = 5;%%to decide how many boxes to consider for each action
    term_f=5;
    term_i_a=4;
    term_f_a=4;
    lambda=[1 0 0];
    high_ordr_exp = 7;
    kthresh=0.25;
    pots_mtrx = ones(num_actions+1)-diag(ones(1,num_actions+1));
    xmld = struct([]);
    
    %%%%%%%%%%%%
    % DISPLAY
    % 0:don't display
    % 1:display the method's results
    % 2:display ground truth results as well
    DISPLAY = 0;
    %%%%%%%%%%%%
    
    %make parfor
    for i=1:num_videos
        ab=videos(i);
        abab = strsplit(ab{1}, '/');
        action = abab{1};
        videoName = abab{2};
        
        flist = dir([img_path '/' action '/' videoName '/*.jpg']);
        num_imgs = length(flist);
        frames = struct([]);
        ts = tic;
        action_paths = cell(num_actions,1);
        imgpathf = [img_path '/'  action '/' videoName];
        xmld(i).videoName = videoName;
        
        %for ground truth
        gtVidInd = getGtVidInd(Videos,videoName);
        gt_tubes = Videos(gtVidInd).tubes;
        
        %%%%%%FOR ssd
        framesssd = readDetections_f([codebase_path , '/datasets/UCF101/ssd-detections/RGB-01-40000/' action '/' videoName '/'],['/home/harkirat/research/16/action/datasets/UCF101/ssd-detections/FLOW-01-75000/' action '/' videoName '/']);
        [action_frames_jpda, Trk, Ff, A_label] = DA_ojla(frames, framesssd, nms_th, num_actions, topk_boxes, imgpathf, action, term_f,gt_tubes,term_f_a,term_i_a, pots_mtrx, lambda, high_ordr_exp, DISPLAY);        
        te = toc(ts);
        
        act_nr = 1;
        min_ac_frames = 20;
        topk = 40; % this is not used in this evaluation
        penalty = 0.5;%for each wrong signal(incorrect action detected)
        l =  strmatch(action,actions,'exact');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% NEW ONLINE LABELLING AND EVALUATION (MULTI-LABEL OJLA)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for a = 1 : size(Trk,2)%nu if tracks
            act_ts = Ff{1,a}(1,1);
            act_te = Ff{1,a}(1,2);
            if act_te<num_imgs
                act_te = act_te - term_f+1;
            end
            bxs=[];
            act_path_scores = [];
            act_path_indices =[];
            boxes=[];
            for t = act_ts:act_te
                bxs = [bxs ; action_frames_jpda(1,t).boxes(a,:)];
                act_path_scores = [act_path_scores ; action_frames_jpda(1,t).scores(a,:)];
                act_path_indices = [act_path_indices ; action_frames_jpda(1,t).actionindices(a)];
            end
            
            bxs = [bxs(:,1:2), bxs(:,3:4)-bxs(:,1:2)+1];
            
            for al=1:num_actions
                label = al;
                act_scores = sort(act_path_scores(:,label),'descend');
                topk_mean = mean(act_scores(1:min(topk,length(act_scores))));

                if (act_te-act_ts+1) >= min_num_frames{l} && topk_mean > kthresh
                    xmld(i).score(act_nr) = double(topk_mean);
                    xmld(i).nr(act_nr) = act_nr;
                    xmld(i).class(act_nr) = label;
                    xmld(i).framenr(act_nr).fnr = act_ts:act_te;

                    xmld(i).boxes(act_nr).bxs = bxs;
                    act_nr = act_nr+1;%no of tracks put in
                end
            end
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% NEW ONLINE-LABELLING (OJLA)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for a = 1 : size(Trk,2)%nu if tracks
%             act_ts = Ff{1,a}(1,1);
%             act_te = Ff{1,a}(1,2);
%             if act_te<num_imgs
%                 act_te = act_te - term_f+1;
%             end
%             bxs=[];
%             path_scores = [];
%             path_indices =[];
%             for t = act_ts:act_te
%                 %you need to take care of a problem of hidden tracks because there's no boxes for it in that frame
%                 path_scores = [path_scores ; action_frames_jpda(1,t).scores(a,:)];
%                 [~,sds] = max(sum(path_scores,1));
%                 action_frames_jpda(1,t).actionindices(a) = sds;
%             end
%             
%         end
%         for a = 1 : size(Trk,2)%nu if tracks
%             act_ts = Ff{1,a}(1,1);
%             act_te = Ff{1,a}(1,2);
%             if act_te<num_imgs
%                 act_te = act_te - term_f+1;
%             end
%             bxs=[];
%             act_path_scores = [];
%             breaks= [act_ts;act_ts;action_frames_jpda(1,act_ts).actionindices(a)];
%             if (act_te-act_ts) > min_num_frames{l}
%                 for t = act_ts:act_te
%                     if action_frames_jpda(1,t).actionindices(a)~= breaks(3,end)
%                         if (breaks(2,end)-breaks(1,end))>=min_num_frames{l} && breaks(3,end)~=25
%     %                         act_scores = sort(act_path_scores,'descend');
%     %                         xmld(i).score(act_nr) = double(mean(act_scores(1:min(topk,length(act_scores)))));
%                             xmld(i).score(act_nr) = double(mean(act_path_scores(:,breaks(3,end))));
%                             xmld(i).nr(act_nr) = act_nr;
%                             xmld(i).class(act_nr) = breaks(3,end);
%                             xmld(i).framenr(act_nr).fnr = breaks(1,end):breaks(2,end);
% 
%                             xmld(i).boxes(act_nr).bxs = [bxs(:,1:2), bxs(:,3:4)-bxs(:,1:2)+1];
%                             act_nr = act_nr+1;%no of tracks put in
% 
%                         end
%                         breaks = horzcat(breaks, [t;t;action_frames_jpda(1,t).actionindices(a)]);
%                         bxs=[];
%                         act_path_scores = [];
%                         bxs = [bxs ; action_frames_jpda(1,t).boxes(a,:)];
%                         act_path_scores = [act_path_scores ; action_frames_jpda(1,t).scores(a,:)];
%                     else
%                         bxs = [bxs ; action_frames_jpda(1,t).boxes(a,:)];
%                         act_path_scores = [act_path_scores ; action_frames_jpda(1,t).scores(a,:)];
%                         breaks(2,end) = t;
%                     end
%                 end
%                 if (breaks(2,end)-breaks(1,end))>=min_num_frames{l} && breaks(3,end)~=25
%                     xmld(i).score(act_nr) = double(mean(act_path_scores(:,breaks(3,end))));
%                     xmld(i).nr(act_nr) = act_nr;
%                     xmld(i).class(act_nr) = breaks(3,end);
%                     xmld(i).framenr(act_nr).fnr = breaks(1,end):breaks(2,end);
% 
%                     xmld(i).boxes(act_nr).bxs = [bxs(:,1:2), bxs(:,3:4)-bxs(:,1:2)+1];
%                     act_nr = act_nr+1;%no of tracks put in
%                 end
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf('done vid=%d\n', i);
    end
    test_list = load([codebase_path , '/datasets/UCF101/ucf101_annot/testlist.mat']);
    testlist=test_list.testlist;
    clear test_list;
    
    %%evaluation
    iou_th=0.2;
    addpath('src/run_evaluation/');
    [mAP,mAIoU, AP] = get_PR_curve_bmvc(Videos, xmld, testlist, actions, iou_th);
    fprintf('mean Average Precision:= %f at spatio-temporal IoU threshold:=%f\n', mAP, iou_th);
    summAP=0;
    for iou_th=0.5:0.05:0.95
        addpath('src/run_evaluation/');
        [mAP,mAIoU, AP] = get_PR_curve_bmvc(Videos, xmld, testlist, actions, iou_th);
        fprintf('mean Average Precision:= %f at spatio-temporal IoU threshold:=%f\n', mAP, iou_th);
        summAP = summAP+mAP;
    end
    fprintf('%f', summAP/10);
    
    %%saving results
    save_path =[codebase_path , '/results'];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    save([save_path '/' 'xmld_' datestr(now,'mm-dd-yyyy HH-MM') '.mat'], 'xmld');
end