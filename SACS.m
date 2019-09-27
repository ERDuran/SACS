%% 
clc
path(pathdef)

birdy = load('chirp');
bird_i_path = audioplayer(birdy.y(1:1500)/5, birdy.Fs/3);
bird_f_path = audioplayer(birdy.y(1500:3000)/5, birdy.Fs/3);

cd ~/SACS
addpath(genpath('~/Dropbox/programs/matlab_functions'))
addpath(genpath('~/Dropbox/programs/ghostscript-9.21'))
addpath(genpath('~/Dropbox/programs/xpdfbin-mac-3.04'))
figures_path = '~/Dropbox/SACS/figures/';
data_path = '~/CARS/';
disp('OK, all set')

set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(groot, 'defaultTextboxshapeInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultPolarAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultTextarrowshapeInterpreter','latex');

