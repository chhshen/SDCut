function [P,lW,bad_indice, sp_maps]=chunk_norris(param)
% function [P,lW,bad_indice]=chunk_norris(param)      
%
% chunk norris :
% He takes pixels and smashes them into superpixels with his bare hands.
%
% lW => nb of superpixels in each image, excluding bad_indice
% P  => projection matrix of the pixels into the superpixel space
% bad_indice => indice of sp without associated subsampled feature points
%
% not used in the article: 
%       support multiscale -- might still contain bugs
%       you can add background images -- might still contain bugs

lW=zeros((param.pic.pics+param.pic.decor)*param.objet.n_lvl_considered,1);


I=[];J=[];


%% OBJ


pos_pt=0;%indice d un pixel, en condiserant toues les images ensemble
ind_pt=0;%indice d un groupe, en bla bla bla

sum_lW=0;

sp_maps = cell(1, param.pic.pics);



for i=1:param.pic.pics

    for j=1:param.objet.n_lvl_considered

        dataBaseSubDir=param.path.obj;
        
        imDir=dir([param.path.im,'superpixel/*_Seg.mat']);
        posDir=dir([dataBaseSubDir,'*_pos']);


        %first open the superpixel groups positions :
        sp_pos=importdata([param.path.im,'superpixel/',imDir(i).name]);
        sp_maps{i} = sp_pos;
        
        %extract group indices :
        sp_ind=unique(sp_pos);

        %position of the param.feature.typeures :
        param.feature.type_pos=importdata([dataBaseSubDir,posDir(i).name]);

        param.feature.type_pos=floor(param.feature.type_pos);
        param.feature.type_pos=param.feature.type_pos(:,2)+size(sp_pos,1)*(param.feature.type_pos(:,1)-1);

        %sparse solution:
        I=[I;(pos_pt+1:pos_pt+size(param.feature.type_pos,1))'];
        J=[J;(sp_pos(param.feature.type_pos)+ind_pt)];

        
        sum_lW=sum_lW+size(sp_ind,1);
        lW((i-1)*param.objet.n_lvl_considered+j)=size(unique(sp_pos(param.feature.type_pos)),1);
        

        ind_pt=ind_pt+size(sp_ind,1);
        pos_pt=pos_pt+size(param.feature.type_pos,1);
        


    end
    
end

%% BACKGROUND -- might still contain bugs

for i=1:param.pic.decor
    for j=1:param.objet.n_lvl_considered
        dataBaseSubDir=param.path.bg;
        
        imDir=dir([param.path.decor,'superpixel/*_Seg.mat']);
        posDir=dir([dataBaseSubDir,'/*_pos']);


        %first open the superpixel groups positions :
        sp_pos=importdata([param.path.decor,'superpixel/',imDir(i).name]);

        %extract group indices :
        sp_ind=unique(sp_pos);
        
        %position of the param.feature.typeures :
        param.feature.type_pos=importdata([dataBaseSubDir,posDir(i).name]);

        param.feature.type_pos=floor(param.feature.type_pos);
        param.feature.type_pos=param.feature.type_pos(:,2)+size(sp_pos,1)*param.feature.type_pos(:,1);

        %sparse solution:
        I=[I;(pos_pt+1:pos_pt+size(param.feature.type_pos,1))'];
        J=[J;(sp_pos(param.feature.type_pos)+ind_pt)];
        sum_lW=sum_lW+size(sp_ind,1);
        
        lW((i-1)*param.objet.n_lvl_considered+j+param.objet.n_lvl_considered*param.pic.pics)=size(unique(sp_pos(param.feature.type_pos)),1);
        

        ind_pt=ind_pt+size(sp_ind,1);
        pos_pt=pos_pt+size(param.feature.type_pos,1);

    end
end

%% CCL


clear param.feature.type_pos sp_pos sp_ind
    
P=sparse(I,J,1,size(I,1),sum_lW);

bad_indice=find(sum(P,1)==0);


P(:,bad_indice)=[];

clear I J

