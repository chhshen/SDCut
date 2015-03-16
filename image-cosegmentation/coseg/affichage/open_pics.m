function pics=open_pics(np,path2,ext,shift,endname)

pics=cell(np,1);

if nargin<5
    endname='';
end

pic_dir=dir([path2,endname,'*',ext]);


i0=1;sh=0;
for i=1:size(pic_dir,1),
    
    if strcmp(pic_dir(i).name(1),'.') || pic_dir(i).isdir
        continue
    end
    
    if sh<shift
        sh=sh+1;
        continue
    end
    
    
    if i0>np
        break
    end
    
    pic_name=[path2,pic_dir(i).name];
    try
        pics{i0}=importdata(pic_name);
    catch
        keyboard
    end
    
    if isstruct(pics{i0})
        try
            pics{i0}=pics{i0}.BW;
        catch
            keyboard
        end
    end
    pics{i0}=double(pics{i0})/255;
    i0=i0+1;

end