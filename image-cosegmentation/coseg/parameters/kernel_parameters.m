% load the parameters for the kernel :
fprintf('# kernels parameters')

% you can either specify directly the constant lambda or the subspace
% dimension df (cf article) :
% ( K + lambda * Id )^-1
param.kernel.df     =100;
param.kernel.lambda =1e-06;


% type of kernel:
param.kernel.type='chi2'; % rbf, linear, chi2 have been implemented...

%cf article
param.kernel.kernel_para=0.1; 

fprintf('\t\t--> done\n')