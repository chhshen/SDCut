function [L]=...
    LaplacianMatrix_pixelwise_between_scales(...
    feat_pos_u,feat_pos_l)


nf_u=size(feat_pos_u,1);% u for upper scale
nf_l=size(feat_pos_l,1);% l for lower scale


if nf_l<nf_u % si la grande scale a plus de features que la petite, on les change
    tmp=feat_pos_u;
    feat_pos_u=feat_pos_l;
    feat_pos_l=tmp;
    clear tmp
end

% Quelque soit l ordre des scales 
% feat_pos_l represente la plus petite des scales 
% donc celle qui a le plus de features

% nbre de connections suivant x et y :

posx_u=unique(feat_pos_u(:,1));
posy_u=unique(feat_pos_u(:,2));
nbx_u=size(posx_u,1);
nby_u=size(posy_u,1);

posx_l=unique(feat_pos_l(:,1));
posy_l=unique(feat_pos_l(:,2));
nbx_l=size(posx_l,1);
nby_l=size(posy_l,1);

%nbre total de connexion :
nb_tot_u=2*(2*nbx_u*nby_u-nby_u-nbx_u);
nb_tot_l=2*(2*nbx_l*nby_l-nby_l-nbx_l);

I=[];J=[];

for i=1:size(feat_pos_u,1)
    
    ligne_i=floor((i-1)/nby_u)+1;
    
    ligne_l=((max(ligne_i*2-1-(ligne_i==1),1):min(ligne_i*2+1+(ligne_i==nbx_u),nbx_l)));
    
    ligne_l=ligne_l*nby_l;
    
    col_i=mod(i,nby_u);
    
    if col_i==0
        col_i=nby_u;
    end
    
    col_l=(max(col_i*2-1-(col_i==1),1):min(col_i*2+1+(col_i==nby_u),nby_l))+1;
    
    try
    ind_l=ones(size(col_l,2),1)*ligne_l+col_l(:)*ones(1,size(ligne_l,2));
    ind_l=ind_l(:)';
    
    
        I=[I;i.*ones(size(ind_l(:)))];
        J=[J;ind_l(:)];

    catch
        keyboard
    end
    
end


%L=sparse(I,J,1);
W=sparse([I;J+size(feat_pos_u,1)],[J+size(feat_pos_u,1);I],1);
clear J I



tmp=1./sum(W).^.5;

I=1:size(tmp(:),1);
D=sparse(I,I,tmp);
Id=sparse(I,I,1);
clear tmp I



L=Id-D*W*D;









