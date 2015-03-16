fprintf('\nLOADING PARAMETERS :\n');
% pic parameters :
pic_parameters;

% object descriptors parameters :
object_parameters;

% optimization parameters :
optimization_parameters;

% paths :
path_parameters;

% kernel parameters :
kernel_parameters;


%obsolete:
param.objet.n_lvl_considered=n_lvl_considered;


%and the number of pics considered :
param.pic.np_considered=param.pic.pics*param.objet.n_lvl_considered;
param.pic.nd_considered=param.pic.decor*param.objet.n_lvl_considered;
fprintf('DONE\n\n');