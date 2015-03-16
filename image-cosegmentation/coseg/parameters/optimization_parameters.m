fprintf('# optimization parameters')



param.optim.sym=0.5;


% weight of the Laplacian matrix :
% ( A + lapWght * L )
if ~exist('lapWght')
	lapWght=0.1;
elseif ischar(lapWght)
    lapWght=str2num(lapWght);
end
param.optim.lapWght=lapWght;



% minimum size of a cluster (in percent / 100 ) :
if ~exist('lambda0')
   lambda0=0.05;
elseif ischar(lambda0) 
    lambda0=str2num(lambda0);
end
param.optim.lambda0=lambda0;


% different mod d optim -- obsolete:
param.optim.mod=1;


% augmented lagrangian param:
param.optim.lambda_init= 1e-03;
param.optim.beta=1.01;


% smoothing =1 <=>  exp(1+log(-alpha*[.]))
param.optim.smoothing=0;
param.optim.alpha=.1;
param.optim.epsilon=0;




fprintf('\t--> done\n')


