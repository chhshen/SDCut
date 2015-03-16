function [R1,R2]=constraints_subroutine_norm3(x,one_pic_bar,one_pic,lambda_1,lambda_2,param)

[R1,R2]=constraints_subroutine(x,one_pic_bar,one_pic,lambda_1,lambda_2,param);
R1=R1.^3;
R2=R2.^3;
    