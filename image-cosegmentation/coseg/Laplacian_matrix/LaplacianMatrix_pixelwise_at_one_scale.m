function [L]=LaplacianMatrix_pixelwise_at_one_scale(feat_pos,Im)
% [L]=LaplacianMatrix_pixelwise_at_one_scale(feat_pos,Im)
%
% compute the Laplacian matrix L for a particular scale
%
% feat_pos: psotion of the pts in the graph
%
% Im: the image




% nb of connections in x and y :

posx=unique(feat_pos(:,1));
posy=unique(feat_pos(:,2));

nbx=size(posx,1);
nby=size(posy,1);

%total nbrof connection :
nb_tot=2*(2*nbx*nby-nby-nbx);


%% evaluate L
I=zeros(nb_tot,1);
J=I;

delta=max(feat_pos(2:end,1)-feat_pos(1:end-1,1));

min_pos=min(feat_pos);
x_min=min_pos(2);y_min=min_pos(1);
max_pos=max(feat_pos);
x_max=max_pos(2); y_max=max_pos(1);

n_pos=size(feat_pos);

subIm=Im(x_min:delta:x_max,y_min:delta:y_max,:);


[i j v] = graphLap(double(subIm),2);

n=numel(subIm)/3;

L=sparse(i+1,j+1,v,n,n);


tmp=1./sum(L).^.5;
I=(1:size(tmp(:)))';
D=sparse(I,I,tmp);
Id=sparse(I,I,1);
clear I tmp


L=Id-D*L*D;




