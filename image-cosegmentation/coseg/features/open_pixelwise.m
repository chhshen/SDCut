function [all_descriptor, all_position,lW,lambda0,feat_size]=open_pixelwise(param)
% function [all_descriptor,all_position,lW,lambda0,feat_size]=open_pixelwise(param)
%
% Open dense descriptors
%
% The descriptor have to be store into a file .mat for each image (for ex:
% knut1_color.mat ) in the directory param.path.obj (for ex:
% ./input/knut/descriptor/color/).
% this file must contain 3 fields: x,y and data.
% x,y gives the positions (they are column vectors with size N)
% data contains the descriptors (each line is a descriptor - therefore it
% is a Nxd matrix)
%
% Not used in the article:
% you can add background pictures for weakly supervised problem

feat_size=0;

subfolder1=num2str(param.objet.patchSize);
dir_names=param.path.obj;
dsifts=dir([dir_names,'*.ma*']);


imageDir=dir([param.path.im,'*.bmp']);


% Then we create a vector of all the descriptors and their location :
lW=zeros(param.pic.np_considered+param.pic.nd_considered,1);% we store the number of descriptors per image and per scale



% evaluate needed memory :
n_desc=zeros(param.objet.n_lvl_considered,2);
n_pos=n_desc;
if param.pic.np_considered>0
    for p=1:param.pic.pics,           
            fsift=importdata([dir_names,dsifts(p).name]);
            n_pos(1)  = n_pos(1)+size(fsift.x,1);
            n_desc(1) = n_desc(1)+size(fsift.data,1);
            if p==1
                n_desc(2) = n_desc(2)+size(fsift.data,2);
            end
            n_pos(2)  = 4;
            clear fsift;
    end
end

all_descriptor = zeros(sum(n_desc(:,1)),n_desc(1,2));
all_position = zeros(sum(n_pos(:,1)),n_pos(1,2));

%%  RETRIEVING OBJECT FEATURES :
fprintf('\n =>loading features:\n')

pos_ind=0;desc_ind=0;

for i=1:param.pic.pics,
        k=1+(i-1)*param.objet.n_lvl_considered;  
        fprintf('\t # %s at scale %i',[dir_names,dsifts(i).name],1)
        fsift=importdata([dir_names,dsifts(i).name]);
        lW(k)=size(fsift.data,1);
        all_position(1+pos_ind:pos_ind+size(fsift.x,1),:)=[floor(fsift.x) , floor(fsift.y) , i*ones(size(fsift.y,1),1), ones(size(fsift.y,1),1)];
        pos_ind=pos_ind+size(fsift.x,1);
        all_descriptor(1+desc_ind:desc_ind+size(fsift.data,1),:)=fsift.data;
        desc_ind=desc_ind+size(fsift.data,1);
        clear fsift
        fprintf('\t\t--> loaded\n')
end
fprintf('   done\n')

desc0=desc_ind;
pos0=pos_ind;



%% BACKGROUND IMAGES (NOT USED IN THE ARTICLE --- might still contain bugs...)


if param.pic.decor~=0

    clear dsifts dir_names
    dir_names=param.path.bg;
    dsifts=dir([dir_names{i},'*ma*']);
    fprintf(' => loading background features:\n')
    for i=1:param.pic.decor,
        k=1+(i-1)*param.objet.n_lvl_considered;
        fprintf('\t # %s at scale %i',[dir_names,dsifts(i).name],1)
        fsift=importdata([dir_names,dsifts(i).name]);
        lW(param.pic.np_considered+k)=size(fsift.data,1);
        all_descriptor=[all_descriptor;fsift.data];
        all_position=[all_position;floor(fsift.x) , floor(fsift.y) , i*ones(size(fsift.y,1),1), ones(size(fsift.y,1),1)];
        clear fsift
        fprintf('\t\t--> loaded\n')
    end
    fprintf('   done\n')
end

%% 

lambda0=floor(lW*param.optim.lambda0+1);