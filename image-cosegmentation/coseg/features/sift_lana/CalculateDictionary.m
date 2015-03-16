function [ ] = CalculateDictionary( imageFileList, dataBaseDir, featureSuffix, dictionarySize, numTextonImages, canSkip , nb , n_lvl )
%function [ ] = CalculateDictionary( imageFileList, dataBaseDir, featureSuffix, dictionarySize, numTextonImages, canSkip )
%
%Create the texton dictionary
%
% First, all of the sift descriptors are loaded for a random set of images. The
% size of this set is determined by numTextonImages. Then k-means is run
% on all the descriptors to find N centers, where N is specified by
% dictionarySize.
%
% imageFileList: cell of file paths
% dataBaseDir: the base directory for the data files that are generated
%  by the algorithm. If this dir is the same as imageBaseDir the files
%  will be generated in the same location as the image files.
% featureSuffix: this is the suffix appended to the image file name to
%  denote the data file that contains the feature textons and coordinates. 
%  Its default value is '_sift.mat'.
% dictionarySize: size of descriptor dictionary (200 has been found to be
%  a good size)
% numTextonImages: number of images to be used to create the histogram
%  bins
% canSkip: if true the calculation will be skipped if the appropriate data 
%  file is found in dataBaseDir. This is very useful if you just want to
%  update some of the data or if you've added new images.

fprintf('Building Dictionary\n\n');

%% parameters

reduce_flag = 1;
ndata_max = 100000;

if(nargin<4)
    dictionarySize = 200
end

if(nargin<5)
    numTextonImages = 50
end

if(nargin<6)
    canSkip = 0
end

if(numTextonImages > size(imageFileList,1) * n_lvl )
    numTextonImages = size(imageFileList,1);
end

outFName = fullfile(dataBaseDir, sprintf('dictionary_%d_%d_%d.mat', dictionarySize, n_lvl,nb ));

if(size(dir(outFName),1)~=0 && canSkip)
    fprintf('Dictionary file %s already exists.\n', outFName);
    return;
end
    

%% load file list and determine indices of training images



inFName = fullfile(dataBaseDir, 'f_order.txt');
if ~isempty(dir(inFName))
    R = load(inFName, '-ascii');
    if(size(R,1)~=size(imageFileList,1))
        R = randperm(size(imageFileList,1));
        sp_make_dir(inFName);
        save(inFName, 'R', '-ascii');
    end
else
    R = randperm(size(imageFileList,1));
    sp_make_dir(inFName);
    save(inFName, 'R', '-ascii');
end

training_indices = R(1:numTextonImages);

%% load all SIFT descriptors

sift_all = [];


for lvl=1:n_lvl
    for f = 1:numTextonImages
        
        imageFName = imageFileList{training_indices(f)};
        [dirN base] = fileparts(imageFName);
        baseFName = fullfile(dirN, base);
        inFName = fullfile(dataBaseDir, sprintf('%d/%d/%s%s', nb*2^(lvl-1),nb ,baseFName, featureSuffix));

        
    
        load(inFName, 'features');
        ndata = size(features.data,1);

        sift_all = [sift_all; features.data];
        fprintf('Loaded %s, %d descriptors, %d so far\n', inFName, ndata, size(sift_all,1));
    end
    
end

fprintf('\nTotal descriptors loaded: %d\n', size(sift_all,1));

ndata = size(sift_all,1);    
if (reduce_flag > 0) & (ndata > ndata_max)
    fprintf('Reducing to %d descriptors\n', ndata_max);
    p = randperm(ndata);
    sift_all = sift_all(p(1:ndata_max),:);
end
        
%% perform clustering
options = foptions;
options(1) = 1; % display
options(2) = 1;
options(3) = 0.1; % precision
options(5) = 1; % initialization
options(14) = 100; % maximum iterations

centers = zeros(dictionarySize, size(sift_all,2));

%% run kmeans
fprintf('\nRunning k-means\n');
dictionary = sp_kmeans(centers, sift_all, options);
    
fprintf('Saving texton dictionary\n');
sp_make_dir(outFName);
save(outFName, 'dictionary');

end
