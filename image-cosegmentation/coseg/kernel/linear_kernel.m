function [xtilde]=linear_kernel(x,df,n,lambda)


xtilde = x - repmat(mean(x,2),1,size(x,2));% COMMENT par A. Joulin : X=Pi_t X



%[u,e] = eig(xtilde * xtilde');
[u,e,v] = svd(xtilde,'econ');
e=e.^2;

df=min(df,floor(size(e,1)/2)+1);

if ~isempty(df);
    lambda = df_to_lambda_eig(df,real(diag(e))/n + 1e-10);
    lambda = 1/lambda;
end

% COMMENT par A. Joulin
% on recupere les valeurs propres de XtX apres centrage superieur a
% n*lambda*0.01 ????

ind = find( diag(e) > n*lambda*1e-2);
q = u(:,ind);
xtilde = q' * xtilde;
% COMMENT par A. Joulin
% C=XX'+n*kappa*Id
C = xtilde * xtilde' + n * lambda * eye(size(xtilde,1));
R = chol(C);%COMMENT par A. Joulin : C=R'R
R = inv(R)';% COMMENT par A. Joulin : R=C^{-1/2}
xtilde = R * xtilde; % COMMENT par A. Joulin : X=C^{-1/2}X => X'X=> X'C^{-1}X

%xtilde=xtilde';%rajout
