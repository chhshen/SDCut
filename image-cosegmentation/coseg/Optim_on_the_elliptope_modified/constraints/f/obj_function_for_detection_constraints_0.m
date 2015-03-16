function f=obj_function_for_detection_constraints_0(x,param,eta)

%  the contribution to the loss function evaluation of the constraints for
%  detection in the case of JBAS / JBAS-I / JBAS-CI
%  



global one_pic_bar % associe a thresh
global one_pic % associe a thresh2

global thresh
global thresh2

global lambda


global old_h
global new_h


f=0;

for i=1:param.pic.pics
    
    lambda_1=thresh(i);% the threshold related to lambda0
    lambda_2=thresh2(i); % **** to n - (k-1) lambda0
    
    if  ~param.optim.smoothing
        [R1_3,R2_3]=constraints_subroutine_norm3(x,one_pic_bar(:,i),one_pic(:,i),lambda_1,lambda_2,param);
    else
        [R1_3,R2_3]=constraints_subroutine_smooth(x,one_pic_bar(:,i),one_pic(:,i),lambda_1,lambda_2,param);
        %R1_3=0.*R1_3;
    end
    
   
 
    
    R1_3=lambdamize(param.optim.norm_bar,lambda,R1_3);
    R2_3=lambdamize(param.optim.norm,lambda,R2_3);

    ff=ones(1,size(x,1))*(R1_3);
    ff=ff+ones(1,size(x,1))*(R2_3);
    
%     if ff~=0
%         keyboard
%     end
    

    f=f+ff;
end

old_h=new_h;
new_h=f;
