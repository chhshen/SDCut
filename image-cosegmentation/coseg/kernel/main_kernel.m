function [xtilde,lW,tab_lambda0,weights]=main_kernel(param)
% function [xtilde,lW,tab_lambda0,weights]=main_kernel(param)
%
% extract features and compute kernel matrix 

% load feature descriptors from saved .mat file
[all_descriptor, trash ,lW,tab_lambda0,trash]=open_pixelwise(param);

clear trash

weights=[sum(lW(1:param.pic.np_considered));sum(lW(param.pic.np_considered+1:end))];

all_descriptor=all_descriptor';

fprintf(' =>evaluating kernel:\n')

[ xtilde , trash,trash ] = compute_kernel(all_descriptor,param.kernel,0,weights);
fprintf('  done')
clear all_descriptor trash