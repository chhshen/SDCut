function [R1,R2]=constraints_subroutine_smooth(x,one_pic_bar,one_pic,lambda_1,lambda_2,param)


[R1,R2]=subroutine_exp(x,one_pic_bar,one_pic,lambda_1,lambda_2,param);

[Rp,Rn]=seperate_signe(R1);
[R1]=subroutine_logExp(Rp,Rn);

[Rp,Rn]=seperate_signe(R2);
[R2]=subroutine_logExp(Rp,Rn);

R1=R1./param.optim.alpha;
R2=R2./param.optim.alpha;
%  R1=(R1./param.optim.alpha).^2;
%  R2=(R2./param.optim.alpha).^2;
 
 