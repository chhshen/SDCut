function plot_groups(...
    param,all_position,z,lW,K,nfig,bad_indice)

color= [ ...
      [1 0 0];...
           [0 1 0];...
       [1 1 0];...
     [0 0.5 .5];... 
     [0 0 1];...
     [0 1 1];...
    [.5 0.5 0];...
    [0.5 0 0];...
    [0.5 0 .5];...
    [1 0 1];...
    [1  1  1];...
    ];

if size(z,2)>1
    z= z*(1:K)';
elseif min(z(:))==0
    z = z+1;
end
    

pos=1;

sumlW=0;
if iscell(param.path.obj)
    paramPathObj = param.path.obj{iObj};
    paramPathIm = param.path.im{iObj};
else
    paramPathObj = param.path.obj;
    paramPathIm = param.path.im;
end
pos_feat_dir= dir([paramPathObj,'*pos']);


for iPic=1:param.pic.pics
        
    pics=open_pics(1,paramPathIm,'',iPic-1);
    pics=pics{1};
    
    
    d=dir([paramPathIm,'superpixel/*_Seg.mat']);
    
    superpix_im=importdata([paramPathIm,'superpixel/',d(iPic).name]);
    
    
    sumlW_future = sumlW + size(unique(superpix_im(:)),1);
    
    superpix_im=superpix_im+sumlW;
    
    
    
    for iBadInd=1:size(bad_indice(:))
        superpix_im(superpix_im==bad_indice(iBadInd))=-1;
    end
    
    allpos=all_position(sum(lW(1:iPic-1))+1: sum(lW(1:iPic)),1 ) + sumlW;
    
    nbox=size(allpos,1);
    
    % image contenant les groupes
    zIm = z(sum(lW(1:iPic-1))+1: sum(lW(1:iPic)) ,:);

    imagek=zeros(size(superpix_im));
    
    for iPix=1:size(allpos,1)
        imagek(superpix_im==allpos(iPix)) = zIm(iPix);
    end
    
    
    pos=pos+nbox;
    
    imagekk = imagek;
    imagek(imagekk(:,1:end-1)~=imagekk(:,2:end,:)) = 0;
    imagekk = imagekk';
    imagek = imagek';
    imagek(imagekk(:,1:end-1)~=imagekk(:,2:end,:)) = 0;
    imagek = imagek';
    imageKcol = reshape(color(imagek+size(color,1)*(imagek==0),:),[size(imagek) 3]);
    
    
    
    
    image_show = 0.4*imageKcol+0.6*(pics.*repmat(imagek~=0,[ 1 1 3]) + repmat(imagek==0,[ 1 1 3]));
    
    
    i0=mod(iPic,16)+16*(mod(iPic,16)==0);
    
    j0=floor((iPic-1)/16);
    figure(nfig+j0)
    subplot (4,4,i0)
    imagesc( image_show );
    set(gca, 'YTick', [])
    set(gca, 'XTick', [])
    hold on
   
    clear filtre_im

    clear pics
    sumlW = sumlW_future;
end



  


end

