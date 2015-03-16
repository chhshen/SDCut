function generate_descriptor_color(param)

% parameters for the dense color descriptors:
gridSpacing     = param.objet.gridSpacing; %1;
patchSize       = param.objet.patchSize;   %1;
assert(gridSpacing == 1);
assert(patchSize == 1);
resolution      = 32;%nb of bins by channel in the histogram

n_lvl           = 1;%for multi-scale purpose...

% path_parameters;

d       = dir([param.path.im,'*',param.objet.type,'*']);
d_size  = size(d,1);

t               = 1;
np_considered   = min(param.pic.pics,d_size);
while t<=np_considered

    if strcmp(d(t).name(1),'.') || d(t).isdir==1
        t = t+1;
        continue
    end
    fprintf('image : \t %s\n',d(t).name)
    
            
    subImageBaseDir = param.path.im; 
    
    imageFileList   = cell(np_considered,1); 
    i0              = 1;
    i               = 1;
    while i0<=np_considered && i<=d_size
        if ~strcmp(d(i).name(1),'.') && d(t).isdir~=1
            imageFileList{i0}   = d(i).name;
            i0                  = i0+1;
        end
        i                       = i+1;
    end


    for i=1:n_lvl
        dataBaseDir_tmp         = param.path.obj;
        GenerateColorDescriptors(imageFileList, subImageBaseDir, dataBaseDir_tmp,gridSpacing/(2^(i-1)), patchSize/(2^(i-1)),resolution);
    end
    
    break
    
end

routine_convert_mat;


end