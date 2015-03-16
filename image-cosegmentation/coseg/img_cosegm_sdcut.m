function img_cosegm_sdcut(in)

type_obj       = in.type_obj;
n_pics         = in.n_pics;

gridSpacing    = in.gridSpacing;
patchSize      = in.patchSize;
feature_type   = in.feature_type;

picMaxSize     = in.picMaxSize;
lapWght        = in.lapWght;
sp_regionSize  = in.sp_regionSize;
sp_regularizer = in.sp_regularizer;

options.fastsdpcut.sigma               = in.sigma;
options.fastsdpcut.lbfgsb_factr        = in.lbfgsb_factr;
options.fastsdpcut.lbfgsb_pgtol        = in.lbfgsb_pgtol;
options.fastsdpcut.lbfgsb_m            = in.lbfgsb_m;
options.fastsdpcut.rounding_method     = in.rounding_method;
options.fastsdpcut.rounding_maxiter    = in.rounding_maxiter;
options.fastsdpcut.rounding_sweepratio = in.rounding_sweepratio;
    
bala_cons_type = in.bala_cons_type;
if bala_cons_type == 2
    lamda  = in.lamda;
else
    lamda0 = in.lamda0;
    lamda1 = in.lamda1;
end
results_name = in.results_name;



parameters;
param.picMaxSize     = picMaxSize;
param.optim.lapWght  = lapWght; %mu in the article 
param.sp_regionSize  = sp_regionSize;
param.sp_regularizer = sp_regularizer;
param.bala_cons_type = bala_cons_type;
if bala_cons_type == 2
    param.lamda  = in.lamda;
else
    param.lamda0 = in.lamda0;
    param.lamda1 = in.lamda1;
end

warning off;

if 1
    if strcmp(feature_type, 'sift')
        generate_descriptor_sift(param);
    else
        generate_descriptor_color(param);
    end

    generate_superpixel(param);
    generate_superpixel_descriptors;
    generate_kernel_laplacian(param);
end

results_dir = [param.path.im, results_name, '/'];
mkdir(results_dir);


%% solve sdp problem
if 1
    data = img_cosegm_sdcut_pack_data(param);
    save([results_dir,'img_cosegm_data.mat'], 'data', 'param', 'options');
else
    load([results_dir,'img_cosegm_data.mat']);
end

if 1
    [~, x_opt, scores, X_opt, u_opt, ~, obj_val1, obj_val2] = img_cosegm_sdcut_solve(data, options.fastsdpcut);
    save([results_dir,'img_cosegm_results.mat'], 'x_opt', 'scores', 'X_opt', 'u_opt', 'obj_val1', 'obj_val2');
else
    load([results_dir,'img_cosegm_results.mat']);
end

%% show results
img_cosegm_plot_result(x_opt,scores,data, param, results_dir);

end

function plot_sp(data, param)

im_dir = param.path.im;
sp_dir = [param.path.im, 'superpixel/'];
sp_suffix = '_Seg.mat';

d             = dir([im_dir,'*',param.objet.type,'*']);
np_considered = min(param.pic.pics , size(d,1));
im_file_list  = cell(np_considered,1);
sp_file_list  = cell(np_considered,1);
for ii=1 : np_considered
    im_file_list{ii} = [im_dir, d(ii).name];
    [~, name, ~] = fileparts(d(ii).name);
    sp_file_list{ii} = [sp_dir, name, sp_suffix];    
end

bad_indicator_full = zeros(1, data.n_sp);
bad_indicator_full(data.bad_indice) = 1;

for ii = 1 : np_considered
    % original img
    subplot(3,np_considered,ii);
    im = imread(im_file_list{ii});
    imshow(im);
    
    % superpixel
    subplot(3,np_considered,ii+np_considered);
    segments.num = data.lW_sp(ii);
    segments.map = data.sp_maps{ii};
    show_segmentation_result(im, segments);
    
    % bad_indice
    subplot(3,np_considered,ii+2*np_considered);
    map = data.sp_maps{ii};
    img = map;
    bad_indicator = bad_indicator_full(data.indice_sp_im{ii});
    for kk = 1 : data.lW_sp(ii)
        if bad_indicator(kk) == 1
            img(map == kk) = -1;
        else
            img(map == kk) = 1;
        end
    end
    imagesc(img); axis image off ; hold on ; 

end

end


function show_segmentation_result(im_rgb, segments)

[sx,sy]=vl_grad(double(segments.map), 'type', 'forward') ;
s = find(sx | sy) ;
imp = im_rgb;
imp([s s+numel(im_rgb(:,:,1)) s+2*numel(im_rgb(:,:,1))]) = 0;

if segments.num == 2
    mask = find(segments.map == 1);
    imp([mask, mask+numel(im_rgb(:,:,1)), mask+2*numel(im_rgb(:,:,1))]) = 0;
end

imagesc(imp) ; axis image off ; hold on ;  

end


















