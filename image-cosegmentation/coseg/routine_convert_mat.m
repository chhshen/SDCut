
obj_name = sprintf('./input/%s/descriptor/' , param.objet.type);



obj_dir=dir(obj_name);
obj_size=size(obj_dir,1);




for o=1:obj_size
    
    if strcmp(obj_dir(o).name(1),'.') || strcmp(obj_dir(o).name,'superpixel')
        continue
    end
    
    feat_name=[obj_name,'/',obj_dir(o).name];
    
    test=dir([feat_name,'/*pos*']);
    
    d=dir([feat_name,'/*.mat']);
    
    for i=1:min(length(d),1000)
        tmp=importdata([feat_name,'/',d(i).name]);
        filename=d(i).name(1:end-4);
        fprintf('processing file : %s \n',filename);
        filename=[feat_name,'/',filename];
        tmp3=tmp.data;
        tmp2=[tmp.x, tmp.y];
        clear tmp
        save('-ASCII',[filename,'_data'],'tmp3');
        save('-ASCII',[filename,'_pos'],'tmp2');
        clear tmp2 tmp3
    end
    
end

