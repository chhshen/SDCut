function [W, aux_data] = calc_affinity_matrix_imgsegm(im, segments, options)

switch options.color_space
    case 'rgb'
        im_adj_mat = im.rgb;
    case 'lab'
        im_adj_mat = im.lab;
    case 'gray'
        im_adj_mat = im.gray;
    case 'hsv'
        im_adj_mat = zeros(size(im.hsv));
        im_adj_mat(:,:,1) = im.hsv(:,:,3);    % v
        im_adj_mat(:,:,2) = im.hsv(:,:,3) .* im.hsv(:,:,2) ...
             .* sin(im.hsv(:,:,1) ./ (2*pi)); % v * s * sin(h)
        im_adj_mat(:,:,3) = im.hsv(:,:,3) .* im.hsv(:,:,2) ...
             .* cos(im.hsv(:,:,1) ./ (2*pi)); % v * s * cos(h)
    otherwise
        error('unknown color space\n')
end

% [h, w, n_ch] = size(im_adj_mat);

% get squared color distance
mean_color    = calc_segments_mean_color(im_adj_mat, segments);
dist_color_sq = calc_pairwise_dist_sq(mean_color);
sigma_color = sqrt(max(dist_color_sq(:))) / 5;
W_color = exp(dist_color_sq / (-1 * sigma_color^2));

% get squared position distance
mean_posit      = calc_segments_mean_pos(segments);
dist_posit_sq   = calc_pairwise_dist_sq(mean_posit);
sigma_posit = sqrt(max(dist_posit_sq(:))) / 10;
W_posit = exp(dist_posit_sq / (-1 * sigma_posit^2));

dist_posit = sqrt(dist_posit_sq);
mink_dist_posit = mink(dist_posit(:), 2*size(dist_posit,1));
mink_dist_posit = mink_dist_posit(size(dist_posit,1)+1:2*size(dist_posit,1));
sigma = mean(mink_dist_posit);


W_posit_threshd = W_posit;
W_posit_threshd(dist_posit_sq > (sigma*4)^2 ) = 0; 

% compute W
W = W_color .* W_posit_threshd;
W = W / max(W(:));
W = (W + W') / 2; 

aux_data.mean_posit = mean_posit;
aux_data.W_color = W_color;
aux_data.W_posit = W_posit;


end




% calculate pairwise distance of a set of colors, 
% a matrix is returned
function dist_sq = calc_pairwise_dist_sq(locations)

[n_points, n_ch] = size(locations);

diff = zeros(n_points, n_points, n_ch);
for ii = 1 : n_ch
    locs_rep1 = repmat(locations(:,ii),  1, n_points);
    locs_rep2 = repmat(locations(:,ii)', n_points, 1);
    diff(:,:,ii) = locs_rep1 - locs_rep2;
end

dist_sq = sum(diff.^2, 3);

end


% calculate mean color of each segment, 
% a vector is returned
function mean_color = calc_segments_mean_color(im_xxx, segments)

[h, w, n_ch] = size(im_xxx);

im_xxx_vec = reshape(im_xxx, [h*w, n_ch]);

mean_color = zeros(segments.num, n_ch);

for ii = 1 : segments.num
    mean_color(ii, :) = sum( im_xxx_vec(segments.map == ii, :), 1);
end

mean_color = mean_color ./ repmat(segments.hist, 1, n_ch);

end

% calculate mean position of each segment
% a vector is returned
function mean_pos = calc_segments_mean_pos(segments)

segments_map = segments.map;

[h, w] = size(segments_map);

mean_pos = zeros(segments.num, 2);

for ii = 1 : w
    for jj = 1 : h
        id = segments_map(jj,ii);
        mean_pos(id, :) = mean_pos(id, :) + [jj,ii];
    end
end

mean_pos = mean_pos ./ repmat(segments.hist, 1, 2);

end




