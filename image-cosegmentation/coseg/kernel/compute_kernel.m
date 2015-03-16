function [ xtilde ,df,H] = compute_kernel(x,kernel_param,somme,weights)
% function [ xtilde ,df,H] = compute_kernel(x,kernel_param,somme,weights)
% 
% compute LINEAR, RBF, HISTORGRAM or CHI-SQUARE kernels
%
% x contains each datapoint as vectors ( dxN matrix)
%
% kernel_param: parameters for the kernel
% somme: useless (warning somme must be equals to 0)

[ d n ] = size(x);

switch kernel_param.type
    case 'linear'
        xtilde =linear_kernel(x,kernel_param.df,n,kernel_param.lambda);
    case 'rbf'
        xtilde = rbf_kernel(x,kernel_param.df,n,kernel_param.kernel_para,kernel_param.lambda);
    case 'kernel_matrix',       
        xtilde = kernel_matrix_kernel(x,kernel_param.df,n);        
    case 'histogram'        
        xtilde = histogram_kernel(x,kernel_param.df,n,kernel_param.lambda);        
    case 'chi2'
        [xtilde,H] = chi2_kernel(x,n,kernel_param.lambda,kernel_param.kernel_para,kernel_param.df,somme,weights);     
end


xtilde=real(xtilde);

df = sum( xtilde(:).^2 );
