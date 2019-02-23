function iou = inters_union(bounds1,bounds2)
    % ------------------------------------------------------------------------
    % [x1 y1 w1 h1] [x2 y2 w2 h2]
    inters = rectint(bounds1,bounds2);
    ar1 = bounds1(:,3).*bounds1(:,4);
    ar2 = bounds2(:,3).*bounds2(:,4);
    union = bsxfun(@plus,ar1,ar2')-inters;
    iou = inters./(union+0.001);
end