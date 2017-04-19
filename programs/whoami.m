%% 
clc
path(pathdef)

Me = input(sprintf([...
    'Which machine is this? \n' ...
    '1: My mac laptop \n' ...
    '2: CCRC linux desktop \n' ...
    '...']));

if Me == 1
    cd ~
    addpath(genpath('Dropbox/SACS_work'))
    addpath(genpath(['/Users/earl/Dropbox/' ...
        'LeeuwinUndercurrent_HonoursProject/' ...
        'matlab/OFAM/ofam_out']))
    outputpath = 'Dropbox/SACS_work/figures/';
    disp('OK, all set')
    
elseif Me == 2
    cd /home/z5100174/Desktop/MATLAB
    addpath(genpath('functions'))
    addpath cars_out
    addpath(genpath('teos10_library'))
    outputpath = 'figures/';
    disp('OK, all set')
    
else
    disp No.
end

