function startup()

    base_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(base_dir, 'functions')));
    addpath(genpath(fullfile(base_dir, 'bin')));
    addpath(genpath(fullfile(base_dir, 'jpda')));
    addpath(genpath(fullfile(base_dir, 'ojla')));
    fprintf('startup done\n');
end