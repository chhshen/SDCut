function img_cosegm(data_name)

addpath(genpath(pwd));

addpath ~/program/PROGLIB/Ncut_9/
addpath ~/program/PROGLIB/lbfgsb.3.0_1.1/lbfgsb3.0_mex1.1
run('~/program/PROGLIB/vlfeat-0.9.14/toolbox/vl_setup');

load('savedState.mat');
stream = RandStream.getGlobalStream;
stream.State = savedState;

run_sdcut(data_name);
run_lowrank(data_name);

end

function run_sdcut(data_name)

if strcmp(data_name, 'car_front')
in.type_obj='car_front';
in.n_pics = 6;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.sigma = 1e-4 / 4017;
in.lbfgsb_factr = 1e7;
in.lbfgsb_pgtol = 1e-5;
in.lbfgsb_m = 200;
in.rounding_method = 'ncut';
in.rounding_maxiter = 10000;
in.rounding_sweepratio = 0.1;
in.bala_cons_type = 2;
in.lamda = 0.9;
in.results_name = 'results_sdcut';
end

if strcmp(data_name, 'car_back')
in.type_obj='car_back';
in.n_pics = 6;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.sigma = 1e-4 / 4012;
in.lbfgsb_factr = 1e7;
in.lbfgsb_pgtol = 1e-5;
in.lbfgsb_m = 200;
in.rounding_method = 'ncut';
in.rounding_maxiter = 10000;
in.rounding_sweepratio = 0.1;
in.bala_cons_type = 2;
in.lamda = 0.9;
in.results_name = 'results_sdcut';
end

if strcmp(data_name, 'face')
in.type_obj='face';
in.n_pics = 10;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.sigma = 1e-4 / 6684;
in.lbfgsb_factr = 1e7;
in.lbfgsb_pgtol = 1e-5;
in.lbfgsb_m = 200;
in.rounding_method = 'ncut';
in.rounding_maxiter = 10000;
in.rounding_sweepratio = 0.1;
in.bala_cons_type = 2;
in.lamda = 0.9;
in.results_name = 'results_sdcut';
end


if strcmp(data_name, 'horse10')
in.type_obj='horse10';
in.n_pics = 10;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.sigma = 1e-4 / 4587;
in.lbfgsb_factr = 1e7;
in.lbfgsb_pgtol = 1e-5;
in.lbfgsb_m = 200;
in.rounding_method = 'ncut';
in.rounding_maxiter = 10000;
in.rounding_sweepratio = 0.1;
in.bala_cons_type = 2;
in.lamda = 0.9;
in.results_name = 'results_sdcut';
end

disp(in.type_obj);
img_cosegm_sdcut(in);

end



function run_lowrank(data_name)


if strcmp(data_name, 'face')
in.type_obj='face';
in.n_pics = 10;
in.gridSpacing = 4; %spacing in the grid
in.patchSize   = 16;%size of a patch
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght=1; %mu in the article -- the only free parameter in our code
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.results_name         = 'results_lowrank';
in.fastsdp_results_name = 'results_sdcut';
end



if strcmp(data_name, 'car_front')
in.type_obj='car_front';
in.n_pics = 6;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.results_name         = 'results_lowrank';
in.fastsdp_results_name = 'results_sdcut';
end

if strcmp(data_name, 'car_back')
in.type_obj='car_back';
in.n_pics = 6;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.results_name         = 'results_lowrank';
in.fastsdp_results_name = 'results_sdcut';
end

if strcmp(data_name, 'horse10')
in.type_obj='horse10';
in.n_pics = 10;
in.gridSpacing = 4;
in.patchSize   = 16;
in.feature_type = 'sift';
in.picMaxSize = 320;
in.lapWght = 1;
in.sp_regionSize = 10;
in.sp_regularizer = 0.5;
in.results_name         = 'results_lowrank';
in.fastsdp_results_name = 'results_sdcut';
end

disp(in.type_obj);
img_cosegm_lowrank(in);

end












