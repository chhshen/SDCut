function generate_superpixel(param)

sp_suffix = '_Seg.mat';

src_dir = param.path.im;
dst_dir = [param.path.im,'superpixel/'];

d             = dir([param.path.im,'*',param.objet.type,'*']);

d_size        = size(d,1);
np_considered = min(param.pic.pics , d_size);
im_list       = cell(np_considered,1);
for i=1:size(d,1)
    im_list{i}=d(i).name;
end

mkdir(dst_dir);

for i=1:np_considered
    im_path = [src_dir, im_list{i}];
    im_ori = single(imread(im_path))./255;
    im = min(max(imresize(im_ori,param.picMaxSize./max(size(im_ori,1),size(im_ori,2))),0),255);
%     [~, NewSeg] = vl_quickseg(im, .7, 2, 15);
    NewSeg = double(vl_slic(im, param.sp_regionSize, param.sp_regularizer) + 1);

    if length(unique(NewSeg)) ~= max(max(NewSeg))
        [NewSeg1, NewSeg1] = histc(NewSeg, unique(NewSeg));
        NewSeg = NewSeg1;
    end
    
    num_sp(i) = max(NewSeg(:));
    [~, im_name, ~] = fileparts(im_list{i}); 
    sp_path = [dst_dir, im_name, sp_suffix];
    save(sp_path, 'NewSeg');
end

sum(num_sp)

end