function [T2,T1]=dual_detection_constraints_subroutine_smooth(x,one_pic_bar,one_pic,lambda_1,lambda_2,param)


[T1,T2]=subroutine_exp(x,one_pic_bar,one_pic,lambda_1,lambda_2,param);

[Rp,Rn]=seperate_signe(T1);
R1=function_sigmExp(Rp,Rn);
T1=subroutine_logExp(Rp,Rn);

[Rp,Rn]=seperate_signe(T2);
R2=function_sigmExp(Rp,Rn);
T2=subroutine_logExp(Rp,Rn);

T1=-2*param.optim.alpha*T1.*R1;
T2=2*param.optim.alpha*T2.*R2;

