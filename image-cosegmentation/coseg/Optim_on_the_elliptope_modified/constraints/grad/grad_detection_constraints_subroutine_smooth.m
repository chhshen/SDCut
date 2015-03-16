function [T1,T2]=grad_detection_constraints_subroutine_smooth(x,...
                                                        one_pic_bar,one_pic,...
                                                        thresh,thresh2,param...
                                                         )
                                                     

[R1,R2]=subroutine_exp(x,one_pic_bar,one_pic,thresh,thresh2,param);

[Rp,Rn]=seperate_signe(R1);
%[T1]=subroutine_logExp(Rp,Rn);
[R1]=function_sigmExp(Rp,Rn);



[Rp,Rn]=seperate_signe(R2);
%[T2]=subroutine_logExp(Rp,Rn);
[R2]=function_sigmExp(Rp,Rn);

% au carre :
% T1 = -2*param.optim.alpha*T1.*R1./(param.optim.alpha*param.optim.alpha);
% T2 = 2*param.optim.alpha*T2.*R2./(param.optim.alpha*param.optim.alpha);

T1 = -param.optim.alpha*R1./(param.optim.alpha);
T2 = param.optim.alpha*R2./(param.optim.alpha);