
randomizer_name=[];
num_ran=size(dir([randomizer_path,'randomizer_*',extension]),1); 

if randomizer~=0 && num_ran==0     
    randomizer_name=[randomizer_path,'randomizer_',num2str(num_ran+1),'_',extension];
else
    tmp=dir([randomizer_path,'randomizer_*',extension]);
    randomizer_name=[randomizer_path,tmp.name];
end





if randomizer~=0
    tmp_dir=dir(randomizer_name);   
    if size(tmp_dir,1)==0
        param.pic.list=randperm(591-33)+33;
        
        param.pic.list=param.pic.list(1:param.pic.pics);
        list=param.pic.list;
        save(randomizer_name,'list');
    else
        param.pic.list=importdata(randomizer_name);
        %param.pic.drift=-1;
    end
end