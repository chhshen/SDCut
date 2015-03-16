function img_cosegm_lowrank(in)

type_obj        = in.type_obj;
n_pics          = in.n_pics;

%param for dense features
gridSpacing     = in.gridSpacing; %spacing in the grid
patchSize       = in.patchSize;%size of a patch
feature_type    = in.feature_type;

picMaxSize      = in.picMaxSize;
lapWght         = in.lapWght; %mu in the article -- the only free parameter in our code
sp_regionSize   = in.sp_regionSize;
sp_regularizer  = in.sp_regularizer;

results_name    = in.results_name;
fastsdp_results_name = in.fastsdp_results_name;

%% parameters

parameters;
param.picMaxSize = picMaxSize;
param.optim.lapWght=lapWght; %mu in the article -- the only free parameter in our code
param.sp_regionSize = sp_regionSize;
param.sp_regularizer = sp_regularizer;

results_dir = [param.path.im, results_name, '/'];
mkdir(results_dir);


%% generate sift descriptors
if 0
generate_descriptor_sift(param);
generate_superpixel(param);
generate_superpixel_descriptors; % WARNOING DOES NOT GENERATE SUPERPIXELS -> only write them to the right format
end


%% FREE PARAMETER TO FIX 


%% Our algorithm
if 1
    y = low_rank_diffrac(param);

    %% Plotting

    global C
    global one_pic

    z=rounding_normalized_cut(y);
    zz=z>0;
    
    

    norm = sqrt(sum(sum((C).^2)));
    C1 = C/norm;
    Y = y*y';
    obj_val1 = Y(:)' * C1(:);
    
    fprintf(1, 'obj_val1=%f\n', obj_val1);
    
    zzz = sign(z);
    obj_val2 = zzz'*C1*zzz;
    fprintf(1, 'obj_val2=%f\n', obj_val2);

    save([results_dir, 'img-cosegm_results_lowrank.mat'], 'C', 'one_pic', 'y', 'z', 'obj_val1', 'obj_val2');
else
    load([results_dir, 'img-cosegm_results_lowrank.mat']);
end
    
    
% plot_detection(zz,param);

scores = z;
x_opt = sign(scores);
x_opt(x_opt == 0) = 1;
fastsdp_results_dir = [param.path.im, fastsdp_results_name, '/'];
fastsdp_data = load([fastsdp_results_dir,'img_cosegm_data.mat']);
data = fastsdp_data.data;

% h = figure;
img_cosegm_plot_result(x_opt,scores,data,param, results_dir);
% print(h, [results_dir, 'img_cosegm_plot_lowrank.eps']);

end






