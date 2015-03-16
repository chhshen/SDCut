function [] = GenerateSiftDescriptors( imageFileList, imageBaseDir, dataBaseDir, maxImageSize, gridSpacing, patchSize, canSkip )
%function [] = GenerateSiftDescriptors( imageFileList, imageBaseDir, dataBaseDir, maxImageSize, gridSpacing, patchSize, canSkip )
%
%Generate the dense grid of sift descriptors for each
% image
%
% imageFileList: cell of file paths
% imageBaseDir: the base directory for the image files
% dataBaseDir: the base directory for the data files that are generated
%  by the algorithm. If this dir is the same as imageBaseDir the files
%  will be generated in the same location as the image files
% maxImageSize: the max image size. If the image is larger it will be
%  resampeled.
% gridSpacing: the spacing for the grid to be used when generating the
%  sift descriptors
% patchSize: the patch size used for generating the sift descriptor
% canSkip: if true the calculation will be skipped if the appropriate data 
%  file is found in dataBaseDir. This is very useful if you just want to
%  update some of the data or if you've added new images.

fprintf('Building Sift Descriptors\n\n');

%% parameters

if(nargin<4)
    maxImageSize = 1000
end

if(nargin<5)
    gridSpacing = 8
end

if(nargin<6)
    patchSize = 16
end

if(nargin<7)
    canSkip = 0
end


for f = 1:size(imageFileList,1)

    %% load image
    imageFName = imageFileList{f};
    
    if isdir(imageFName)
        continue
    end
    
    [dirN base] = fileparts(imageFName);
    baseFName = [dirN filesep base];
    outFName = fullfile(dataBaseDir, sprintf('%s_sift.mat', baseFName));
    imageFName = fullfile(imageBaseDir, imageFName);
    
    if(size(dir(outFName),1)~=0 && canSkip)
        fprintf('Skipping %s\n', imageFName);
        continue;
    end
    
    I = sp_load_image(imageFName);

    [hgt wid] = size(I);
    if min(hgt,wid) > maxImageSize
        I = imresize(I, maxImageSize/min(hgt,wid), 'bicubic');
        fprintf('Loaded %s: original size %d x %d, resizing to %d x %d\n', ...
            imageFName, wid, hgt, size(I,2), size(I,1));
        [hgt wid] = size(I);
    end

    %% make grid (coordinates of upper left patch corners)
    remX = mod(wid-patchSize,gridSpacing);
    offsetX = floor(remX/2)+1;
    remY = mod(hgt-patchSize,gridSpacing);
    offsetY = floor(remY/2)+1;
    
    [gridX,gridY] = meshgrid(offsetX:gridSpacing:wid-patchSize+1, offsetY:gridSpacing:hgt-patchSize+1);

    fprintf('Processing %s: wid %d, hgt %d, grid size: %d x %d, %d patches\n', ...
             imageFName, wid, hgt, size(gridX,2), size(gridX,1), numel(gridX));

    %% find SIFT descriptors
    siftArr = sp_find_sift_grid(I, gridX, gridY, patchSize, 0.8);
    siftArr = sp_normalize_sift(siftArr);
    
    features.data = siftArr;
    features.x = gridX(:) + patchSize/2 - 0.5;
    features.y = gridY(:) + patchSize/2 - 0.5;
    features.wid = wid;
    features.hgt = hgt;

    sp_make_dir(outFName);

    save(outFName, 'features');

end % for

end % function
