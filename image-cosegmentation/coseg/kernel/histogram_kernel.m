function [xtilde]=histogram_kernel(x,df,n,lambda)
x=x';
H=zeros(size( x * x' ) );

for i=1:size(H,1)
    for j=1:size(H,2)
        H(i,j)=sum(min(x(i,:)'.*x(i,:)',x(j,:)'.*x(j,:)'));
        %H(i,j)=sum(min(x(i,:)',x(j,:)'));
    end
end



%G=cholinc(H,'0.1');
[V,D]=eig(H);

G=V*sqrt(D);

G=G';

G = G - repmat(mean(G,2),1,size(G,2));% COMMENT par A. Joulin : X=Pi_t X

C = G * G' +  n * lambda * eye(size(G,1));

try
    R = chol(C);
catch
    % augment lambda if not positive enough and error is produced
    R = chol(C + norm(G,'fro')^2 * 1e-10 * eye(size(G,1)));
    warning('lambda was augmented');
end

 R = inv(R)';

xtilde = R * G;


