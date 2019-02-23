function [boxes, scores] = read_dt_files(dt_boxes_file)
    st = load(dt_boxes_file);
    boxes = st.boxes;
    scores =  st.scores;
    clear st;
end