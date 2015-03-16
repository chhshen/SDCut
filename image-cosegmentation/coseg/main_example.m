clear all 
cd /home/peng/program/FastSDPCut/image-cosegmentation/coseg
addpath(genpath(pwd));

if 0
    type_obj='coke';
    n_pics = 2;

    %param for dense features
    % gridSpacing = 4; %spacing in the grid
    % patchSize   = 16;%size of a patch
    gridSpacing = 1; %spacing in the grid
    patchSize   = 1;%size of a patch
    feature_type = 'color';

    parameters;
    param.picMaxSize = 128;
    param.optim.lapWght=0.001; %mu in the article -- the only free parameter in our code
    param.sp_regionSize = 10;
    param.sp_regularizer = 0.5;
end

if 1
    type_obj='face';
    n_pics = 10;

    %param for dense features
    gridSpacing = 4; %spacing in the grid
    patchSize   = 16;%size of a patch
    feature_type = 'sift';

    picMaxSize = 320;
    lapWght=1; %mu in the article -- the only free parameter in our code
    sp_regionSize = 10;
    sp_regularizer = 0.5;
end

%% parameters

parameters;
param.picMaxSize = picMaxSize;
param.optim.lapWght=lapWght; %mu in the article -- the only free parameter in our code
param.sp_regionSize = sp_regionSize;
param.sp_regularizer = sp_regularizer;


%% generate sift descriptors
if 1

generate_descriptor_sift(param);

generate_superpixel(param);

generate_superpixel_descriptors; % WARNOING DOES NOT GENERATE SUPERPIXELS -> only write them to the right format
end


%% FREE PARAMETER TO FIX 


%% Our algorithm
y = low_rank_diffrac(param);

%% Plotting

global C
global one_pic

z=rounding_normalized_cut(y);
zz=z>0;

save('lowrank-data.mat', 'C', 'one_pic', 'y', 'z');

plot_detection(zz,param);







