function imgsegm_biased(im_file, options)
%% read image
global im;
im.rgb = im2single(imread(im_file));
im.lab = im2single(vl_xyz2lab(vl_rgb2xyz(im.rgb)));
im.lab = im.lab ./ max(abs(im.lab(:)));
im.gray = rgb2gray(im.rgb);
im.hsv = rgb2hsv(imread(im_file));
im.h = size(im.rgb, 1);
im.w = size(im.rgb, 2);

%% over segmentation
switch options.over_segm.color_space
    case 'rgb'
        im_over_segm = im.rgb;
    case 'lab'
        im_over_segm = im.lab;
    case 'gray'
        im_over_segm = im.gray;
    otherwise
        error('unknown color space: %s\n', options.color_space);
end
% SEGMENTS = VL_SLIC(IM, REGIONSIZE, REGULARIZER) 
global over_segments;
over_segments.map  = vl_slic(im_over_segm, options.over_segm.region_size, 0.1);
over_segments.map  = over_segments.map + 1;
over_segments.num  = max(over_segments.map(:)); %length(unique(over_segments_map));
over_segments.hist = single(histc(over_segments.map(:), 1:over_segments.num));


%% obtain constrains
pos_in = options.pos_in;
neg_in = options.neg_in;

h = size(over_segments.map, 1);

pos_pixels_idxs = (pos_in(:,1)-1) * h + pos_in(:,2);
bias_group_pos = unique(over_segments.map(pos_pixels_idxs));

if options.num_neg_in > 0
    neg_pixels_idxs = (neg_in(:,1)-1) * h + neg_in(:,2);
    bias_group_neg = unique(over_segments.map(neg_pixels_idxs));
else
    bias_group_neg = [];
end

options.same_pairs = [];
options.diff_pairs = [];
options.bias_kappa_same = options.bias_kappa_same;
options.bias_kappa_diff = options.bias_kappa_diff;
options.bias_group_pos = bias_group_pos;
options.bias_group_neg = bias_group_neg;

options.sigma = options.sigma / double(over_segments.num);

%% compute affinity matrix for over-segments
global aux_data;
[W, aux_data] = calc_affinity_matrix_imgsegm(im, over_segments, options.adj_mat);


%% binary cut on graph of over-segments
for ii = 1 : length(options.graph_cut_methods)
    switch options.graph_cut_methods{ii}
        case 'sdcut'
            [~, x, x_real] = solve_sdcut(sparse(W), options);
        case 'bncut'
            [~, x, x_real] = solve_bncut(sparse(W), options);
        otherwise 
            error('unknown options.cut_method: %s\n', options.cut_method);
    end

    segments.map  = ismember(over_segments.map, find(x==1));
    segments.map  = segments.map + 1;
    segments.num  = 2;
    segments.hist = single(histc(segments.map, 1:2));


    show_segmentation_result_and_score(im.rgb, segments, options, ii);
end


end

%%
function show_segmentation_result_and_score(im_rgb, segments, options, method_ii)

figure;
imp = im_rgb;
if segments.num == 2
    map = segments.map-1;
    pos_label = 0;
    for ii = 1 : size(options.pos_in,1)
        pos_label = pos_label + map(options.pos_in(ii,2), options.pos_in(ii,1));
    end
    if pos_label/size(options.pos_in,1) < 0.5
        map = 1 - map;
    end
    
    map = bwareaopen(map, 500);
    map = 1-bwareaopen(1-map, 500);
        
    mask = find(map ~= 1);
    imp([mask, mask+numel(im_rgb(:,:,1)), mask+2*numel(im_rgb(:,:,1))]) = 0;
end
imagesc(imp) ; axis image off ; hold on ;  
% saveas(gcf, [options.rpath, options.im_name, '/', options.im_name, '_', options.graph_cut_methods{method_ii}, '_segm_result.fig']);
saveas(gcf, [options.rpath, options.im_name, '/', options.im_name, '_', options.graph_cut_methods{method_ii}, '_segm_result.eps'], 'psc2');


end











