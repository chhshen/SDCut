function subgrad=grad_detection_constraints_0(x,param,eta)

%  the contribution to the gradient evaluation of the constraints for
%  detection in the case of JBAS / JBAS-I / JBAS-CI
% 


global one_pic_bar % associe a thresh
global one_pic

global thresh
global thresh2

global lambda


sum_x=sum_x_over_V(x,one_pic(:,1));
sum_x_bar=sum_x_over_V(x,one_pic_bar(:,1));
if ~param.optim.smoothing
    [T1,T2]=grad_detection_constraints_subroutine(x,one_pic_bar(:,1),one_pic(:,1),thresh(1),thresh2(1),param);
else
    [T1,T2]=grad_detection_constraints_subroutine_smooth(x,one_pic_bar(:,1),one_pic(:,1),thresh(1),thresh2(1),param);
end
T1=lambdamize(param.optim.norm_bar,lambda,T1);
T2=lambdamize(param.optim.norm,lambda,T2);

subgrad=(T2*sum_x+one_pic(:,1)*(T2'*x));
subgrad=subgrad+(T1*sum_x_bar+one_pic_bar(:,1)*(T1'*x));

for i=2:param.pic.pics 
    
    sum_x=sum_x_over_V(x,one_pic(:,i));
    sum_x_bar=sum_x_over_V(x,one_pic_bar(:,i));
    
    if ~param.optim.smoothing
        [T1,T2]=grad_detection_constraints_subroutine(x,one_pic_bar(:,i),one_pic(:,i),thresh(i),thresh2(i),param);
        T1=3*T1;
        T2=3*T2;
    else
        [T1,T2]=grad_detection_constraints_subroutine_smooth(x,one_pic_bar(:,i),one_pic(:,i),thresh(i),thresh2(i),param); 
    end
    
    T1=lambdamize(param.optim.norm_bar,lambda,T1);
    T2=lambdamize(param.optim.norm,lambda,T2);
    
    subgrad=subgrad+(T2*sum_x+one_pic(:,i)*(T2'*x));
    subgrad=subgrad+(T1*sum_x_bar+one_pic_bar(:,i)*(T1'*x));
    
end