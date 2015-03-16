
featureSuffix   = '_sift.mat'; 


path_parameters;

param.path.obj  = [param.path.im,'descriptor/sift/'];


d               = dir([param.path.im,'*',param.objet.type,'*']);
d_size          = size(d,1);




np_considered   = min(param.pic.pics , d_size);

imageFileList   = cell(np_considered,1);

for i=1:size(d,1)
    imageFileList{i}=d(i).name;
end

maxImageSize    = 1000;


GenerateSiftDescriptors( imageFileList, param.path.im, param.path.obj, maxImageSize, gridSpacing,patchSize, 0 );
clear param.path.feat_tmp


routine_convert_mat;
