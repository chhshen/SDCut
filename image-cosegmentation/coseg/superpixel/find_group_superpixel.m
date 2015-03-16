function find_group_superpixel(type_obj,imageFileList,superpixFileList,dataBaseDir,imageBaseDir, patchSize,npmax)

superpixDir=[imageBaseDir,type_obj,'/desc_superpixel/'];

[dirPath trash] = fileparts(superpixDir);
if(isdir(dirPath)==0)
    mkdir(dirPath);
end

np                  = min(size(imageFileList,1),npmax);%nb pics


for ff = 1:np
    f=ff;
    fprintf('image %i\n',ff)
   
    %on cherche les superpixels :
    tmp2 = importdata([imageBaseDir,superpixFileList{ff}]);
    [ind,pos_ind]=unique(tmp2);
    %pour chaque boite on va chercher les sift qui tombe dedans :
    pos=zeros(size(ind,1),3);


    for i=1:size(ind,1)

        pos(i,1)=tmp2(pos_ind(i));
        pos(i,2)=sum(sum(tmp2==ind(i)));
        pos(i,3)=f;


    end


    outFileName=[superpixDir,type_obj,sprintf('_%03d',f)];
    outFileName=[outFileName,'_pos'];
    save(outFileName,'pos');
    fprintf(['saved in ',outFileName,'\n'])

end

