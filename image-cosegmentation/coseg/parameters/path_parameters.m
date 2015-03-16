% parameters for paths :
fprintf('# paths') 


data_path = '.';

bg_folder='bg';

%where to save things...
save_path=[data_path,'/output/'];

%where to import things...
import_path=[data_path,'/input/'];

% images location
param.path.im = sprintf('%s%s/',import_path,param.objet.type);



if size(dir(param.path.im),1)==0
    fprintf('\n\nERROR: input images not found\n\n');
    keyboard
end

% if background images are used:
param.path.decor= sprintf('%s%s/',import_path,bg_folder);



% obj descriptor location :
param.path.obj = [param.path.im,  'descriptor/', feature_type, '/'];
% background pics location :
param.path.bg = [param.path.decor,'descriptor/', feature_type, '/'];

%% saving files :

%extension des fichier


subfields=fieldnames(param);


save_path=[save_path,param.objet.type,'/'];
save_path_file=['y_nbIm',num2str(param.pic.pics)];

param.path.save_file=save_path_file;
sp_make_dir( save_path );

param.path.save=save_path;





fprintf('\t\t\t\t--> done\n')