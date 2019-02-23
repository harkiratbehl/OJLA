function [action_frames_jpda, Trk, Ff, A_label] = DA_ojla(frames, framesssd, nms_th, num_actions, topk_boxes ,imgpathf,action, term_f, gt_tubes,term_f_a,term_i_a, pots_mtrx, lambda, high_ordr_exp, DISPLAY)
        
    a=2;
    flist = dir([imgpathf '/*.jpg']);%Boxes have coordinates as x1,y1,x2,y2
    num_imgs = length(flist);
    a = a+1; % to skip background boxes and scores
    action_frames = struct([]);
    action_frames_jpda = struct([]);
    actions = {'Basketball', 'BasketballDunk', 'Biking', 'CliffDiving', 'CricketBowling', 'Diving', 'Fencing',...
    'FloorGymnastics', 'GolfSwing', 'HorseRiding', 'IceDancing', 'LongJump', 'PoleVault', 'RopeClimbing', 'SalsaSpin', ...
    'SkateBoarding', 'Skiing', 'Skijet', 'SoccerJuggling', 'Surfing', 'TennisSwing', 'TrampolineJumping', 'VolleyballSpiking', 'WalkingWithDog'};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PARAMETERS
    cost0=10;
    mbest=2;
    topk_boxes=3;
    topk_actions=10;
    min_init_score=0.1;
    term_f=5;
    min_score_to_assoc=0.08;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_gt_tubes = length(gt_tubes);
        
    for f = 1:length(framesssd)   
        fname = flist(f).name;
        frimg = [imgpathf  '/' fname];
        ab = imread(frimg);
        if DISPLAY ~= 0
            imshow(ab)
            hold on
        end 
        %%------ SPATIAL-----------
        boxes_spatial = framesssd(f).catboxes;
        scores_spatial = framesssd(f).catscores;
        [sortedX,sortingIndicess] = sort(scores_spatial,2,'descend');
        score = sortedX(:,1:topk_actions);
        boxes_spatial_mixed = boxes_spatial;
        scores_spatial_mixed = zeros(size(scores_spatial,1),1);
        for l=1:size(boxes_spatial,1)
            [scores_spatial_mixed(l,1) ~] = max(scores_spatial(l,1:end-1));
        end
        boxes_spatial = [boxes_spatial_mixed, scores_spatial_mixed];
        boxes_spatial_before = boxes_spatial_mixed;
        scores_spatial_before = scores_spatial;
        
        pick_nms_spatial = nms(boxes_spatial, nms_th);    
        boxes_spatial_all = boxes_spatial_before(pick_nms_spatial,:);
        scores_spatial_all = scores_spatial_before(pick_nms_spatial,:);
        
        if length(pick_nms_spatial)>topk_boxes
            pick_nms_spatial = pick_nms_spatial(1:topk_boxes); 
        end
        boxes_spatial = boxes_spatial_before(pick_nms_spatial, :); 
        scores_spatial = scores_spatial_before(pick_nms_spatial,:);

        
        %%%%%%%%%%%%%%%%%%%%
        keep = [];
        for b=1:size(boxes_spatial_all,1)
            [a, in]=max(scores_spatial_all(b,:));
            if a>min_score_to_assoc && in~=25
                Xx = double(boxes_spatial_all(b,1));
                Yy = double(boxes_spatial_all(b,2));
                Ww = double(boxes_spatial_all(b,3)-Xx);
                Hh = double(boxes_spatial_all(b,4)-Yy);
%                 rectangle('Position',[Xx,Yy,Ww,Hh],'EdgeColor',[0,0,0],'LineWidth',2)
%                 text(Xx+Ww/2,Yy+Hh/3,num2str(scores_spatial_all(b,1)),'Color',[1,1,1],'FontSize',15) %action score
    %                 text(Xx+Wnum_gt_tubes = length(gt_tubes); w/2,Yy+Hh/2,actions(a),'Color',[1,1,1],'FontSize',15) %action score
                keep = [keep ; b];
            end
        end
        boxes_spatial_all = boxes_spatial_all(keep,:);
        scores_spatial_all = scores_spatial_all(keep,:);
        
        boxes_spatial_all = boxes_spatial;
        scores_spatial_all = scores_spatial;
        
        
        %% OJLA
        if (f==1)
            boxes_spatial_all = boxes_spatial;
            scores_spatial_all = scores_spatial;            
            action_frames_jpda(f).boxes = boxes_spatial_all;
            action_frames_jpda(f).scores = scores_spatial_all;
            [sortX,sortIndices] = max(scores_spatial_all,[],2);
            action_frames_jpda(f).actionindices = sortIndices;
            
            % Frme has the no of frames in total
            MTch_Tg=cell(1,length(framesssd));%make a cell of matrix of size no. of frames(one matrix for each frame)
            F_Pr=cell(1,length(framesssd));
            Ntg=size(action_frames_jpda(f).boxes,1);%no. of bounding boxes in 1st frame(initial no. of tracks)
            Trk=cell(1,Ntg);%for each track
            Exist_trg=1:Ntg;%it stores the labels for all the existing tracks
            Oc_cnt=zeros(1,Ntg);
            %%%
            Action_in_cnt=zeros(1,Ntg);
            Action_oc_cnt=zeros(1,Ntg);
            A_label=cell(1,Ntg);%final actual label
            A_label_ins=cell(1,Ntg);%instantaneous best label label
            %%
            Ff=cell(1,Ntg);%stores start and end frame index(last in which it's seen) for each track
            for i=Exist_trg%i is the track label
                Trk{i}=[i;0;1];%put these into track cell for each track...[3] stores the frame no. [1] stores it's bounding box no. in that frame
                Ff{1,i}=[1 1];
                %%%
                A_label{i}=[0];
                A_label_ins{i} = action_frames_jpda(f).actionindices(i);
                %%%
            end
            color=.25+.75*rand(Ntg,3);
            colorbbnum=.25+.75*rand(1,3);
%             Gate = (boxes_spatial(1,3)-boxes_spatial(1,1))/2;
            Gate=100;
    
        else
            
            D_si2=size(boxes_spatial_all,1);%no. of bounding boxes in kth frame
            MTch_Tg{f}=cell(max(Exist_trg),1);%create a matrix for this frame which is itself a cell of size of tracks
            Mes_Tar=false(D_si2,max(Exist_trg));
            inside_gate=false(D_si2,1);
            Mtch_indx = zeros(D_si2,max(Exist_trg));%stores the index of matching 
            for i=Exist_trg%over all the existing tracks
                MTch_Tg{f}{i}.Costs(1)=cost0;
                MTch_Tg{f}{i}.Meas_edge=[0;1];%for the no detection case
                if Trk{i}(1,end)~=0%if it was matched to something in the previous frame
                    xy=[(action_frames_jpda(f-1).boxes(i,1)+action_frames_jpda(f-1).boxes(i,3))/2 (action_frames_jpda(f-1).boxes(i,2)+action_frames_jpda(f-1).boxes(i,4))/2];%puts the centre box
                    bounds1=[action_frames_jpda(f-1).boxes(i,1) action_frames_jpda(f-1).boxes(i,2) (-action_frames_jpda(f-1).boxes(i,1)+action_frames_jpda(f-1).boxes(i,3)+1) (-action_frames_jpda(f-1).boxes(i,2)+action_frames_jpda(f-1).boxes(i,4)+1)];
                    lnz=f-1;
                else
                    lnz=find(Trk{i}(1,:)~=0); % use this for stored features
                    lnz = Trk{i}(3,lnz(end));
                    xy=[(action_frames_jpda(lnz).boxes(i,1)+action_frames_jpda(lnz).boxes(i,3))/2 (action_frames_jpda(lnz).boxes(i,2)+action_frames_jpda(lnz).boxes(i,4))/2];
                    bounds1=[action_frames_jpda(lnz).boxes(i,1) action_frames_jpda(lnz).boxes(i,2) (-action_frames_jpda(lnz).boxes(i,1)+action_frames_jpda(lnz).boxes(i,3)+1) (-action_frames_jpda(lnz).boxes(i,2)+action_frames_jpda(lnz).boxes(i,4)+1)];
                end
                
                %%%%% displaying for better debugging
%                 rectangle('Position',bounds1,'EdgeColor',[1 1 1],'LineWidth',2)
%                 text(double(bounds1(1)+bounds1(3)/2),double(bounds1(2)+bounds1(4)/3),num2str(i),'Color',[1 1 1],'FontSize',18)
                %%%%
                
                % 1. computing the higher order term matrix
                nz=find(Trk{i}(1,:)~=0);
                nz = Trk{i}(3,nz);
                joint_t = zeros(1,num_actions+1);
                if length(nz)>=high_ordr_exp
                    nz = nz(end-high_ordr_exp+1:end);
                    for wp=1:length(nz)
                        joint_t = joint_t + action_frames_jpda(nz(wp)).scores(i,:);
                    end
                    for l=1:num_actions+1
                        joint_t(l) = joint_t(l) - pots_mtrx(action_frames_jpda(lnz).actionindices(i),l);%potts term
                    end
                else
                    for wp=1:length(nz)
                        joint_t = joint_t + action_frames_jpda(nz(wp)).scores(i,:);
                    end%haven't decided action yet if we haven't found it till now
                end
                 
                 
                Gate = (action_frames_jpda(f-1).boxes(i,3)-action_frames_jpda(f-1).boxes(i,1))/2;
                % Gatr size has to increase with mis-detections
                
                for j=1:D_si2%over all the detected bounding boxes in this frame
                    if pdist2(xy,[(boxes_spatial_all(j,1)+boxes_spatial_all(j,3))/2 (boxes_spatial_all(j,2)+boxes_spatial_all(j,4))/2])<=Gate
                        
                        %%DECIDING WHICH (higher order terms and potts potential)
                        joint_term = joint_t + scores_spatial_all(j,:);
                        [joint_term, indx] = max(joint_term);
                        
                        Mtch_indx(i,j) = indx;
                        
                        %BINARY POTENTIALS
                        bounds2=[boxes_spatial_all(j,1) boxes_spatial_all(j,2) (-boxes_spatial_all(j,1)+boxes_spatial_all(j,3)+1) (-boxes_spatial_all(j,2)+boxes_spatial_all(j,4)+1)];
                        costiou = inters_union(bounds1,bounds2);
                        disc = pdist2(xy,[(boxes_spatial_all(j,1)+boxes_spatial_all(j,3))/2 (boxes_spatial_all(j,2)+boxes_spatial_all(j,4))/2]);
                        dotp = 0;
                        
                        cos = double(1/(joint_term+((length(nz)+1)/2)*lambda(1)*costiou+lambda(2)*disc+lambda(3)*dotp));
                        
                        MTch_Tg{f}{i}.Costs = [MTch_Tg{f}{i}.Costs;cos];
                        MTch_Tg{f}{i}.Meas_edge=[MTch_Tg{f}{i}.Meas_edge [j;1]];%all the matching's that can exist
                        Mes_Tar(j,i)=true;
                        inside_gate(j,1)=true;
                    end
                end

                MTch_Tg{f}{i}.Hypo=MTch_Tg{f}{i}.Meas_edge(1,:)';%target i has been matched(computed costs) with these detections in kth frame
                MTch_Tg{f}{i}.Prob=exp(-MTch_Tg{f}{i}.Costs);%converting costs into probabilities
                MTch_Tg{f}{i}.A_Eq_Const=sparse(ones(1,size(MTch_Tg{f}{i}.Meas_edge,2)));
                MTch_Tg{f}{i}.b_Eq_Const=1;
            end
            ixxd=(~cellfun(@isempty,MTch_Tg{f}));%puts 1 for the existing tracks only
            F_Pr{f}=cell(1,size(MTch_Tg{f},1));%cell of size max(Exist tracks)
            Mes_Tar2=Mes_Tar(:,ixxd);[Umt, Vmt]=size(Mes_Tar2);
            Mes_Tar3=[false(Vmt,Vmt+Umt);Mes_Tar2 false(Umt,Umt)];
            if mbest>0
                F_Pr{f}(ixxd) =Approx_Multiscan_JPDA_Probabilities(Mes_Tar3,MTch_Tg{f}(ixxd),mbest);%pass only the existing tracks data
            elseif mbest==0
                F_Pr{f}(ixxd)=cellfun(@(x) x.Prob/sum(x.Prob),MTch_Tg{f}(ixxd),'UniformOutput', false);
            else
                error('mbest value should be equal or bigger than zero')
            end
%             hold off; drawnow;
%             pause(0.1)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% FOR TRACK UPDATION(1st markov chain)
            assigned_meas=[];
            for i=Exist_trg
                [Vlu,Ixx]=max(F_Pr{1,f}{i});
                Ix=MTch_Tg{f}{i}.Hypo(Ixx);
                Trk{i}=[Trk{i} [Ix;Vlu;f]];
                if Ix~=0
                    assigned_meas=[assigned_meas Ix]; %#ok<AGROW>
                    % Exist_trg_new=[Exist_trg_new i]; %#ok<AGROW>
                    Oc_cnt(i)=0;
                    action_frames_jpda(f).boxes(i,:) = boxes_spatial_all(Ix,:);
                    action_frames_jpda(f).scores(i,:) = scores_spatial_all(Ix,:);
                    action_frames_jpda(f).actionindices(i) = Mtch_indx(i,Ix);
                else
                    Oc_cnt(i)=Oc_cnt(i)+1;
                    action_frames_jpda(f).boxes(i,:) = action_frames_jpda(f-1).boxes(i,:);
                    action_frames_jpda(f).scores(i,:) = action_frames_jpda(f-1).scores(i,:);
                    action_frames_jpda(f).actionindices(i) = action_frames_jpda(f-1).actionindices(i);
                end
            end
            
                     
%% TRACK INITIATION
            %initiate only those which are out of gate of anyone and in
            %top 2 for some action
            n_trgt=setdiff(1:D_si2,assigned_meas);
            trc_siz=size(Trk,2);
            new=0;
            for i=1:length(n_trgt)
                ij=n_trgt(i);
                [X,Indic] = max(scores_spatial_all(ij,:));
                %conditions
                if inside_gate(ij,1)~=true && X>min_init_score && Indic~=25
                    new=new+1;
                    Trk{trc_siz+new}=[ij;0;f];
                    color=[color;.25+.75*rand(1,3)];
                    Oc_cnt(trc_siz+new)=0;
                    action_frames_jpda(f).boxes(trc_siz+new,:) = boxes_spatial_all(ij,:);
                    action_frames_jpda(f).scores(trc_siz+new,:) = scores_spatial_all(ij,:);
                    action_frames_jpda(f).actionindices(trc_siz+new) = Indic;
                    Ff{1,trc_siz+new}=[f f];
                    A_label{trc_siz+new}=[0];
                    A_label_ins{trc_siz+new} = Indic;
                    Action_in_cnt(trc_siz+new) = [0];
                    Action_oc_cnt(trc_siz+new) = [0];
                end
            end
            Exist_trg_new=1:size(Trk,2);

%% TRACK TERMINATION
            Exist_trg=Exist_trg_new(Oc_cnt<term_f);
            for i=Exist_trg
                Ff{1,i}(1,2)=f;
            end
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%% THIS PART DISPLAYS ACTIVE TRACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if DISPLAY~=0
            for i=Exist_trg
                kt=find(Trk{i}(3,:)==f);
                if Trk{i}(1,kt)~=0
                    Xx = double(action_frames_jpda(f).boxes(i,1));
                    Yy = double(action_frames_jpda(f).boxes(i,2));
                    Ww = double(action_frames_jpda(f).boxes(i,3)-Xx+1);
                    Hh = double(action_frames_jpda(f).boxes(i,4)-Yy+1);
                    rectangle('Position',[Xx,Yy,Ww,Hh],'EdgeColor',color(i,:),'LineWidth',2)
    %                 text(Xx+Ww/2,Yy+Hh/4,num2str(A_label{i}(end)),'Color',color(i,:),'FontSize',15) %action score
                    text(Xx+Ww/2,Yy+Hh/2,num2str(action_frames_jpda(f).actionindices(i)),'Color',color(i,:),'FontSize',15) %action score
                    text(Xx+Ww/2,Yy+Hh/10,num2str(i),'Color',color(i,:),'FontSize',20)
                end
            end
            text(size(ab,2)/2,size(ab,1)-10,num2str(f),'FontSize',25)

            %displaying ground truth
            if DISPLAY==2
                for gtind = 1:num_gt_tubes
                    gt_bb = gt_tubes(gtind).boxes;
                    gt_fnr = gt_tubes(gtind).sf:gt_tubes(gtind).ef;
                    if f>=gt_tubes(gtind).sf && f<=gt_tubes(gtind).ef
                        Xx = double(gt_bb(f-gt_fnr(1)+1,1));
                        Yy = double(gt_bb(f-gt_fnr(1)+1,2));
                        Ww = double(gt_bb(f-gt_fnr(1)+1,3));
                        Hh = double(gt_bb(f-gt_fnr(1)+1,4));
                        rectangle('Position',[Xx+1,Yy+1,Ww+3,Hh+3],'EdgeColor',[0 0 0],'LineWidth',2)
                    end
                end
            end
            hold off; drawnow;
            pause(0.01)
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    clear boxes scores pick_nms pick_softmax;
end