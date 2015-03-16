function y = low_rank_diffrac(param)
% function low_rank_diffrac(param)
% MAIN FUNCTION :
% Use the low rank optimization algorithm by Journ\'ee et al.
% it is based on the diffrac framework (Bach and Harchoui '07)
% 
% Its purpose is to detect a common object in a set of pictures



% we are optimizing a max-cut function over the elliptope :
fun_set             = @functions_elliptope;
fun_obj             = @functions_max_cut;
fun_constraints     = @functions_constraints_0;


%% COMPUTE KERNEL

global one_pic
global one_pic_bar
global C


fprintf('COMPUTE KERNEL...')
tic
[xtilde,lW,param.optim.tab_lambda0,weights]     = main_kernel(param);
t   = toc;
fprintf('\n Elapsed time for the computation of the kernel is %f seconds\nDONE\n\n',t)

[ d n ] = size(xtilde);


%% Laplacian Matrix

fprintf('COMPUTE LAPLACIAN MATRIX L...\n')
tic
L   = LaplacianMatrix_pixelwise(param,lW,n);
t   = toc;
fprintf(' Elapsed time for the computation of L is %f seconds\nDONE\n\n',t)

%% Chunks into SUPERPIXELS

[one_pic]           = find_indices_image(lW,ones(sum(lW),1));
one_pic_bar         = one_pic;

%% Compute the Discriminative clustering matrix A

fprintf('COMPUTE DISCRIMINATIVE CLUSTERING MATRIX A...\n')
tic
lW_px               = lW;
[P,lW,bad_indice]   = chunk_norris(param);
C                   = compute_A(param,xtilde,L,lW_px,lW,weights,P);
clear xtilde L


one_pic             = P'*one_pic;
one_pic_bar         = P'*one_pic_bar;
clear P

%% POSITIVE CONSTRAINTS (needed for background -- might contain bugs) see
%% Bach and Harchoui (2007) for that matter.

if param.pic.nd_considered ~= 0
    positive_constraints                       = sum(lW(1:param.pic.np_considered))+1;
    size_pos_const                             = sum(lW(param.pic.np_considered+1:end));
    one_pic(positive_constraints,:)            = sum(one_pic(positive_constraints:end,:));
    one_pic                                    = one_pic(1:positive_constraints,:);
    one_pic(:,end-param.pic.nd_considered+1)   = sum(one_pic(:,end-param.pic.nd_considered+1:end),2);
    one_pic                                    = one_pic(:,1:end-param.pic.nd_considered+1);
    one_pic_bar                                = one_pic;
    
    I  = (1:sum(lW))';
    J  = [(1:sum(lW(1:param.pic.np_considered)))';...
          (sum(lW(1:param.pic.np_considered))+1)*ones(sum(lW(param.pic.np_considered+1:end)),1)];
    PC = sparse(I,J,1);
    clear I J
    
    C  = PC'*C*PC;
    clear PC
    
    lW      = lW(1:param.pic.np_considered+(param.pic.nd_considered~=0));
    lW(end) = lW(end)*(param.pic.nd_considered==0)+(param.pic.nd_considered~=0);
    
    param.pic.tab_lambda0       = param.pic.tab_lambda0(1:end-param.pic.nd_considered+(nd_considered~=0));
    param.pic.tab_lambda0(end)  = param.pic.tab_lambda0(end)*(param.pic.nd_considered==0);
end


p   = sum(lW);

trC = trace(C);
C   = C/trC;
t   = toc;

fprintf(' Elapsed time for the computation of A is %f seconds\nDONE\n\n',t)


%% CONSTRAINTS

global thresh
global thresh2

% transform for the elliptope :
n                       = size(C,1);
thresh                  = (2*param.optim.tab_lambda0-sum(one_pic)');
thresh2                 = (sum(one_pic)'-2*param.optim.tab_lambda0);

param.optim.norm        = normalization(param.pic.pics,sum(one_pic(:)));
param.optim.norm_bar    = normalization(param.pic.pics,sum(one_pic_bar(:)));

%% Low-rank optimization :

fprintf('LOW RANK OPTIMIZATION...\n')
tic
% INITIALIZATION:
r0      = 2;
y0      = feval(fun_set,'retraction',randn(p,r0),zeros(p,r0),param);
% OPTIMIZATION:
[y,f]   = low_rank_optim(fun_set,fun_obj,param,y0,'TR',fun_constraints);
t       = toc;
fprintf(' Elapsed time for the optimization is %f seconds\n',t)
fprintf('DONE\n\n')
