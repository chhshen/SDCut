function GenerateColorDescriptors(imageFileList, imageBaseDir, dataBaseDir,gridSpacing, patchSize,resolution)
%function [] = GenerateSiftDescriptors( imageFileList, imageBaseDir, dataBaseDir, maxImageSize, gridSpacing, patchSize, canSkip )
%
% Generate the dense grid of color histogram for each
% image
% the size of an histogram is patchSize^2

fprintf('Building Color Descriptors\n\n');


for f = 1:size(imageFileList,1)
    
    if size(imageFileList{f},1)==0
        break
    end
    
    %% load image
    imageFName = imageFileList{f};
    [dirN base] = fileparts(imageFName);
    baseFName = [dirN filesep base];
    
    outFName = fullfile(dataBaseDir, sprintf('%s_color.mat', baseFName));
    imageFName = fullfile(imageBaseDir, imageFName);
    
    I = double(importdata(imageFName));
    I = I/256;
    
    [hgt wid trash] = size(I);
    
    
    %% make grid (coordinates of upper left patch corners)
    remX = mod(wid-patchSize,gridSpacing);
    offsetX = floor(remX/2)+1;
    remY = mod(hgt-patchSize,gridSpacing);
    offsetY = floor(remY/2)+1;
    
    [gridX,gridY] = meshgrid(offsetX:gridSpacing:wid-patchSize+1, offsetY:gridSpacing:hgt-patchSize+1);
    %[gridX,gridY] = meshgrid(offsetX:gridSpacing:wid-16+1, offsetY:gridSpacing:hgt-16+1);
    
    fprintf('Processing %s: wid %d, hgt %d, grid size: %d x %d, %d patches\n', ...
             imageFName, wid, hgt, size(gridX,2), size(gridX,1), numel(gridX));
         

    features=sp_find_color_grid(gridX,gridY,I, patchSize,resolution);
    
    sp_make_dir(outFName);

    save(outFName, 'features');

    
end

