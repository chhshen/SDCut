function L=LaplacianMatrix_pixelwise(param,lW,n)
% function L=LaplacianMatrix_pixelwise(param,lW,n)

% EVALUTATE THE LAPLACIAN MATRIX
% when we use dense descriptors on a gride.

% NOT IN THE ARTICLE -- can still contain bugs:
% 
% It can deal with multiple scales (quad tree).
% The output is the combination of the Laplacian matrices at different
% scales.

np_considered=param.pic.np_considered;

I=[];
J=[];
V=[];


subfolder1='/color';

L_size=0;

imdir=dir([param.path.im,param.objet.type,'*']);

for np=1:np_considered/param.objet.n_lvl_considered     
    
    fprintf(' # image %i...',np)
    old_L_size=L_size;    
    
    for i=1:param.objet.n_lvl_considered        

        dataBaseSubDir=param.path.obj;        
%         posDir=dir([dataBaseSubDir,'/*_pos']);
        posDir=dir([dataBaseSubDir,'*_pos']);

        % to find the need size for the Laplacian matrix L for an image:
%         param.feature.type_pos=importdata([dataBaseSubDir,'/',posDir(np).name]);
        
param.feature.type_pos=importdata([dataBaseSubDir,posDir(np).name]);
        param.feature.type_pos=floor(param.feature.type_pos);
        
        Im=importdata([param.path.im,imdir(np).name]);
        Im=double(Im);
        
        % Laplacian for a scale :
        L=LaplacianMatrix_pixelwise_at_one_scale(param.feature.type_pos,Im);
        
        [I_tmp,J_tmp,V_tmp]=find(L);
   
        I=[I;I_tmp+L_size];
        J=[J;J_tmp+L_size];
        V=[V;V_tmp];
        
        L_size=L_size+size(L,1);   
        clear L
    end    
    
    
    % Laplacian between scales --might still contain bugs:    
    for i=1:param.objet.n_lvl_considered-1       
        for j=i+1:param.objet.n_lvl_considered
            subfolder1=num2str(param.objet.patchSize*j);
            dataBaseSubDir=[param.path.obj,'/',subfolder1];
            posDir=dir([dataBaseSubDir,'/*_pos']);
            param.feature.type_pos_2=floor(importdata([dataBaseSubDir,'/',posDir(1).name]));
            
            L2=LaplacianMatrix_pixelwise_between_scales(param.feature.type_pos,param.feature.type_pos_2);
            
            [I_tmp,J_tmp,V_tmp]=find(L2);
            I=[I;I_tmp+old_L_size];
            J=[J;J_tmp+old_L_size];
            V=[V;V_tmp];
            
            clear L2
        end
    end
    
    fprintf(' done\n')
end


% background -- might still contain bugs:
for np=1:param.pic.nd_considered/param.objet.n_lvl_considered 
    old_L_size=L_size;    
    for i=1:param.objet.n_lvl_considered        
        if param.feature.type==0
            subfolder1=['/',num2str(param.objet.patchSize*i)];
        end
        dataBaseSubDir=[param.path.bg,subfolder1];        
        posDir=dir([dataBaseSubDir,'/*_pos']);
        param.feature.type_pos=importdata([dataBaseSubDir,'/',posDir(np).name]);
        param.feature.type_pos=floor(param.feature.type_pos);        
        % Laplacian inside a scale :
        L=LaplacianMatrix_pixelwise_at_one_scale(param.feature.type_pos);        
        [I_tmp,J_tmp,V_tmp]=find(L);
        I=[I;I_tmp+L_size];
        J=[J;J_tmp+L_size];
        V=[V;V_tmp];
        L_size=L_size+size(L,1);
        clear L
    end
    
    % Laplacian between scales :
    for j=i+1:param.objet.n_lvl_considered
        subfolder1=num2str(param.objet.patchSize*j);
        dataBaseSubDir=[param.path.bg,'/',subfolder1];
        posDir=dir([dataBaseSubDir,'/*_pos']);
        param.feature.type_pos_2=floor(importdata([dataBaseSubDir,'/',posDir(1).name]));
        L2=LaplacianMatrix_pixelwise_between_scales(param.feature.type_pos,param.feature.type_pos_2);
        [I_tmp,J_tmp,V_tmp]=find(L2);
        I=[I;I_tmp+old_L_size];
        J=[J;J_tmp+old_L_size];
        V=[V;V_tmp];
        
        clear L2
    end
end

L=sparse(I,J,V);
clear I J V

if size(L,1)~=n
    L=kron(eye(param.pic.pics),L);
    L=sparse(L);
end


