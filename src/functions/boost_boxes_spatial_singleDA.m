function [bs, s_s, bf ,s_f] = boost_boxes_spatial_singleDA(bs, bf, s_s, s_f) % bs - boxes_spatial bf-boxes_flow
    nb = size(bs,1); % num boxes
    iou_thr = 0.3;

    box_spatial = [bs(:,1:2) bs(:,3:4)-bs(:,1:2)+1];
    box_flow =    [bf(:,1:2) bf(:,3:4)-bf(:,1:2)+1];


    for i=1:nb
        ovlp = inters_union(box_spatial(i,:), box_flow); % ovlp has 1x5 or 5x1 dim
        [movlp, mind] = max(ovlp);
        if movlp>=iou_thr
            for j=1:size(s_s,2)
                s_s(i,j) = (s_s(i,j) + s_f(mind,j)*movlp)/2;
            end
        end
    end
end