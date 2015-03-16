function img_cosegm_plot_result(x_opt, scores, data, param, results_dir)


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

u_full = zeros(1, data.n_sp);
u_full(data.good_indice) = scores;

x_full = zeros(1, data.n_sp);
x_full(data.good_indice) = x_opt;


if 1
      
for ii = 1 : np_considered

    
    % scores
    h = figure;
    map = data.sp_maps{ii};
    if length(unique(map)) ~= max(max(map))
        [map1, map1] = histc(map, unique(map));
        map = map1;
    end
    img = map;
    u = u_full(data.indice_sp_im{ii});
    for kk = 1 : data.lW_sp(ii)
        img(map == kk) = u(kk);
    end
    imagesc(-1 * img); axis image off; hold on; 

    
    filename = [results_dir, param.objet.type, '_score_fastsdp_', num2str(ii)];

    saveas(h, filename, 'fig');
    saveas(h, [filename, '.eps'], 'psc2');
end
end


if 0
    
im_sizes = [250 197
            225 150
            198 220
            309 195
            275 184
            304 279
            220 183
            220 178
            220 181
            220 187];    
    
for ii = 1 : np_considered
    im_size = im_sizes(ii,:);
    
    % scores
    h = figure;
    map = data.sp_maps{ii};
    if length(unique(map)) ~= max(max(map))
        [map1, map1] = histc(map, unique(map));
        map = map1;
    end
    img = map;
    u = u_full(data.indice_sp_im{ii});
    for kk = 1 : data.lW_sp(ii)
        img(map == kk) = u(kk);
    end
    imagesc(-1 * img(1:im_size(2),1:im_size(1))); axis image off; hold on; 

    
    filename = [results_dir, param.objet.type, '_score_fastsdp_', num2str(ii)];

    saveas(h, filename, 'fig');
    saveas(h, [filename, '.eps'], 'psc2');
end
end

if 0
    h = figure;
    for ii = 1 : np_considered
        % original img
        subplot(4,np_considered,ii);
        im = imread(im_file_list{ii});
        image(im); axis image off; hold on;

        % superpixel
        subplot(4,np_considered,ii+np_considered);
        image(label2rgb(data.sp_maps{ii}, rand(max(data.sp_maps{ii}(:)),3), [0 0 0], 'shuffle')); axis image off; hold on;

        % scores
        subplot(4,np_considered,ii+2*np_considered);
        map = data.sp_maps{ii};
        if length(unique(map)) ~= max(max(map))
            [map1, map1] = histc(map, unique(map));
            map = map1;
        end
        img = map;
        u = u_full(data.indice_sp_im{ii});
        for kk = 1 : data.lW_sp(ii)
            img(map == kk) = u(kk);
        end
        imagesc(img); axis image off; hold on; 


        % x
        subplot(4,np_considered,ii+3*np_considered);
        map = data.sp_maps{ii};
        if length(unique(map)) ~= max(max(map))
            [map1, map1] = histc(map, unique(map));
            map = map1;
        end
        img = map;
        x = x_full(data.indice_sp_im{ii});
        for kk = 1 : data.lW_sp(ii)
            img(map == kk) = x(kk);
        end
        imagesc(img); axis image off ; hold on ; 

        %         map = double(over_segments.map);
    %         for kk = 1 : n
    %             map(map == kk) = u(kk);           
    %         end
    %         imagesc(map); axis image off ; hold on ; 
    %         title('order');



    %     idx_start = x_opt_idx;
    %     idx_end   = x_opt_idx + max(NewSeg(:))-1;
    %     sp_idx_pos = find(x_opt(idx_start:idx_end) == 1);
    %     segments.num = 2;
    %     segments.map = ones(size(NewSeg));
    %     segments.map(ismember(NewSeg,sp_idx_pos)) = 2;
    %     show_segmentation_result(im, segments);   
    %     x_opt_idx = idx_end + 1;
    end
    print(h, [results_dir, 'img_cosegm_plot.eps']);
end

end