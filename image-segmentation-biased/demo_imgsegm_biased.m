function demo_imgsegm_biased()

%% add paths
%% Please change to the correct path on your machine
run('~/program/PROGLIB/vlfeat-0.9.14/toolbox/vl_setup');
addpath ~/program/PROGLIB/lbfgsb.3.0_1.1/lbfgsb3.0_mex1.1
addpath ~/program/PROGLIB/MinMaxSelection/MinMaxSelection;

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);

clear all;
close all;
demo_imgsegm_biased_i('55073');
clear all;
close all;
demo_imgsegm_biased_i('69020');
clear all;
close all;
demo_imgsegm_biased_i('113044');
clear all;
close all;
demo_imgsegm_biased_i('161045');
clear all;
close all;
demo_imgsegm_biased_i('207056');

end

function demo_imgsegm_biased_i(im_name)


%% input paras
options.im_name = im_name;
options.rpath = './results/';

im_file = ['./images/', options.im_name, '.jpg'];


options.bias_kappa_same = 0.75; % larger  than
options.bias_kappa_diff = 0.75; % smaller than

%
options.balance_cons_type = 'non-equpartition';
options.a1 = 0.4;
options.a2 = 0;

%
options.over_segm.color_space = 'lab';
options.over_segm.region_size = 15;

%
options.adj_mat.calc_method = 'Ncut-PAMI';
options.adj_mat.color_space = 'lab';

%
options.graph_cut_methods = {'bncut', 'sdcut'};


options.sigma = 1e-2;
options.lbfgsb_factr = 1e7;
options.lbfgsb_pgtol = 1e-5;
options.lbfgsb_m = 200;
options.rounding_method = 'srh';

options.cvx_solver = '';
options.cvx_precision  = 'default';

options.biased_ncut_nvec   = 26;


%% Set/Load biased points
options.num_pos_in = 10;
options.num_neg_in = 10;
    
num_pos_in = options.num_pos_in;
num_neg_in = options.num_neg_in;

if 0
    % To set the bias points
    I = imread(im_file);
    [nr,nc,nb] = size(I);


    figure(1); clf;
    imshow(I); hold on;

    fprintf(1, 'click %i pos points\n', num_pos_in);
    
    pos_in = [];
    for ii = 1 : num_pos_in
        pos_in_single = round(ginput(1));
        plot(pos_in_single(1), pos_in_single(2), 'r+', 'MarkerSize', 8, 'LineWidth', 3); hold on;
        pos_in(ii,:) = pos_in_single;
    end
    

    if options.num_neg_in > 0
        fprintf(1, 'click %i neg points\n', num_neg_in);
        
        neg_in = [];
        for ii = 1 : num_neg_in
            neg_in_single = round(ginput(1));
            plot(neg_in_single(1), neg_in_single(2), 'bx', 'MarkerSize', 8, 'LineWidth', 3);
            neg_in(ii,:) = neg_in_single;
        end
    else
        neg_in = [];
    end
    save([options.rpath, options.im_name, '/', options.im_name, '_biased_points.mat'], 'pos_in', 'neg_in');
    
    keyboard
    
else
    % To load the bias points
    load([options.rpath, options.im_name, '/', options.im_name, '_biased_points.mat']);

    I = imread(im_file);
    figure(1); clf;
    image(I); axis image off ; hold on;
    
    plot(pos_in(:,1), pos_in(:,2), 'r+', 'MarkerSize', 8, 'LineWidth', 3);
    if options.num_neg_in > 0
        plot(neg_in(:,1), neg_in(:,2), 'bx',  'MarkerSize', 8, 'LineWidth', 3);
    else
        neg_in = [];
    end
end

saveas(gcf, [options.rpath, options.im_name, '/', options.im_name, '_img_with_bias.eps'], 'psc2');

options.pos_in = pos_in;
options.neg_in = neg_in;

%%
imgsegm_biased(im_file, options);

end
