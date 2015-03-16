function [xtilde,H]=chi2_kernel(x,n,lambda,chi2_lambda,df,somme,weights)



H=[];

if ~isempty(df)
    mmax=min(1000,4*df);
else
    mmax=1000;
end

fprintf('\t C code running...')

if somme==0
    [G,P,m] = icd_chi2(x,chi2_lambda,n*1e-8,min(n,mmax)); 
else
    [G,P,m] = icd_chi2_somme(x,chi2_lambda,n*1e-8,min(n,mmax),somme);
end


fprintf('done\n')

[trash,Pi]=sort(P);
clear P
G=G(Pi,1:m);
clear Pi


w_pos=0;
for i=1:size(weights,1)
    G(w_pos+1:w_pos+weights(i),:)=G(w_pos+1:w_pos+weights(i),:)./sqrt(weights(i));
end

G=G';

%H=G'*G;

G = G - repmat(mean(G,2),1,size(G,2));

[u,e,v] = svd(G,'econ');
e=e.^2;
ind =find( real(diag(e))/n > 1e-10 );
if df >= length(ind),
    % df is too big
    lambda = 1e-10;
else
    lambda = df_to_lambda_eig(df,real(diag(e(ind,ind)))/n);
    lambda = 1/lambda;
end
ind = find( real(diag(e)) > n*lambda*1e-2);
q = u(:,ind);
G = q' * G;

C = G * G' +  lambda * eye(size(G,1));
%C = G * G' +  n * lambda * eye(size(G,1));

try
    C = chol(C);
    
catch
    % augment lambda if not positive enough and error is produced
    warning('lambda was augmented');
    keyboard
end
C = inv(C)'; % r * r matrix

xtilde = C * G; % r * n matrix
clear C G


