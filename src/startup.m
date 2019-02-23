function startup()

    base_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(base_dir, 'utils')));
    addpath(genpath(fullfile(base_dir, 'functions')));
    addpath(genpath(fullfile(base_dir, 'bin')));
    addpath(genpath(fullfile(base_dir, 'jpda')));

    mkdir_if_missing(fullfile(base_dir, 'external'));

    caffe_path = fullfile(base_dir, 'external', 'caffe', 'matlab');
    if exist(caffe_path, 'dir') == 0
        error('matcaffe is missing from external/caffe/matlab; See README.md');
    end
    addpath(genpath(caffe_path));


    fprintf('fast_rcnn startup done\n');
end