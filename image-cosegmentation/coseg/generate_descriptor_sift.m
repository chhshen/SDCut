function generate_descriptor_sift(param)

% featureSuffix   = '_sift.mat'; 


param_path_obj  = [param.path.im,'descriptor/sift/'];


d               = dir([param.path.im,'*',param.objet.type,'*']);
d_size          = size(d,1);

np_considered   = min(param.pic.pics, d_size);

imageFileList   = cell(np_considered,1);

for i=1:np_considered
    imageFileList{i}=d(i).name;
end

maxImageSize    = 1000;


GenerateSiftDescriptors( imageFileList, param.path.im, param_path_obj, maxImageSize, ...
                         param.objet.gridSpacing, param.objet.patchSize, 0 );

routine_convert_mat;

end