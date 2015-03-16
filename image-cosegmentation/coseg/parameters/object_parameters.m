% set parameters specific to an object 
% ( dense sift specifications ...)


% for dense sift/color hist
%         -->  gridSpacing (4 usually)
%         -->  patchSize (16 for sift)


% for multi scale descriptors 
%         --> n_lvl 

fprintf('# object specifications')



n_lvl=1;
n_lvl_considered=max(1,n_lvl);

%% DEFINE PARAM.OBJET

param.objet.type=type_obj;

param.objet.patchSize=patchSize;

param.objet.gridSpacing=gridSpacing;

param.objet.n_lvl=n_lvl;


% if superpixels are used
if ~exist('superpix')
    superpix=1;
elseif ischar(superpix)
    superpix=str2num(superpix);
end
param.feature.superpix=superpix;


fprintf('\t\t--> done\n')