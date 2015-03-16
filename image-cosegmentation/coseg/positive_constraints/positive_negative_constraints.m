function [xtilde,D,PC,L,one_pic,one_pic_bar] ...
        =positive_negative_constraints(xtilde,...
                                        D,n, ...
                                        L,...
                                        positive_constraints, ...
                                        one_pic,one_pic_bar,n_decor...
                                        )


PC=[];

p=size(D,1);

 if ~isempty(positive_constraints)
    % POSITIVE CONSTRAINTS
    % transform pairs of matching points to chunks
    
    np=size(positive_constraints,1);
    
    D(p-np)=sum(D(p-np:end));
    D=D(1:p-np);
    
    one_pic(p-np,:)=sum(one_pic(p-np:end,:));
    one_pic=one_pic(1:p-np,:);
    one_pic(:,end-n_decor+1)=sum(one_pic(:,end-n_decor+1:end),2);
    one_pic=one_pic(:,1:end-n_decor+1);
    
    one_pic_bar(p-np,:)=sum(one_pic_bar(p-np:end,:));
    one_pic_bar=one_pic_bar(1:p-np,:);    
    one_pic_bar=one_pic_bar(:,1:end-n_decor+1);
    one_pic_bar(end,end)=0;

    I=[(1:p-np-1)';(p-np:p)'];
    J=[(1:p-np-1)'; (p-np).*ones(np+1,1)];
    PC=sparse(I,J,1);


    xtilde = xtilde * PC;
    try
    L       =PC'*L*PC;
    catch
    tmp=zeros(positive_constraints(1,2));
    tmp(1:positive_constraints(1,2)-1,1:positive_constraints(1,2)-1)=L(1:positive_constraints(1,2)-1,1:positive_constraints(1,2)-1);
    tmp(end)=positive_constraints(end,1)-positive_constraints(end,2)+1;
    L=tmp;
    clear tmp
    end

 end
 

