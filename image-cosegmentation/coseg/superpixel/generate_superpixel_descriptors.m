
imageBaseDir=['./input/'];
% imageBaseDir=['/home/peng/program/FastSDPCut/image-cosegmentation/coseg/input/'];

dataBaseDir=sprintf('./input/%s/descriptor/',param.objet.type);
% dataBaseDir=sprintf('/home/peng/program/FastSDPCut/image-cosegmentation/coseg/input/%s/descriptor/',param.objet.type);

d=dir([imageBaseDir,type_obj,'/superpixel/*_Seg.mat']);
dd=dir([imageBaseDir,type_obj,'/',type_obj,'*']);
dd=dd(1:length(d));


superpixFileList=cell(size(d));
imageFileList=superpixFileList;
ii=1;



for i=1:size(dd,1)
    
    if strcmp(dd(i).name(1),'.') || isdir([imageBaseDir,type_obj,'/',dd(i).name])
        continue
    end
    imageFileList{ii}=[type_obj,'/',dd(i).name];
    
    try
        superpixFileList{ii}=[type_obj,'/superpixel/',d(ii).name];
    catch
        keyboard
    end
    ii=ii+1;    
end

find_group_superpixel(type_obj,imageFileList,superpixFileList,dataBaseDir ,imageBaseDir,patchSize,size(d,1));


