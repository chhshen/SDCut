function [all_descriptor, all_position,lW,tab_lambda0,bad_ind]=open_superpix(param)

pos_obj=dir([param.path.im,'desc_superpixel/',param.objet.type,'*pos*']);



tab_lambda0=zeros(param.pic.np_considered+param.pic.decor,1);

all_descriptor=[];
all_position=[];
lW=[];

sumTmp=0;
bad_ind=[];


for i=1:param.pic.pics
    for k=1:param.objet.n_lvl_considered

        
        pos_obj=dir([param.path.im,'desc_superpixel/',param.objet.type,sprintf('_%03i',i),'*pos*']);

       
        tmp2=importdata([param.path.im,'desc_superpixel/',pos_obj(1).name]);
       
        tmp2=[tmp2, (k).*ones(size(tmp2,1),1) ];
        
        all_position=[all_position;tmp2];
        tab_lambda0(i)=max(floor(param.optim.lambda0*size(tmp2,1)),1);
        sumTmp=sumTmp+size(tmp2,1);
        lW=[lW;size(tmp2,1)];

    end
end
