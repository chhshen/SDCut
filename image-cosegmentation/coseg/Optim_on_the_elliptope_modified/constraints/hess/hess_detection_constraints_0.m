function subhess=hess_detection_constraints_0(x,param,eta)

%  the contribution to the hessian evaluation of the constraints for
%  detection in the case of JBAS / JBAS-I / JBAS-CI
% 

global npics


global one_pic_bar
global one_pic

global thresh
global thresh2

global lambda



subhess=zeros(size(x));

for i=1:param.pic.pics
    lambda_1=thresh(i);
    lambda_2=thresh2(i);

    
    if ~param.optim.smoothing
        tmp=hess_subroutine_0(x,eta,one_pic_bar(:,i),lambda_1,param.optim.epsilon,-1);
        tmp=lambdamize(param.optim.norm_bar,lambda,tmp);
        subhess=subhess+tmp;
        tmp=hess_subroutine_0(x,eta,one_pic(:,i),lambda_2,param.optim.epsilon,1);
        tmp=lambdamize(param.optim.norm_bar,lambda,tmp);       
        subhess=subhess+tmp;
    else
        tmp=hess_subroutine_smooth(x,eta,one_pic_bar(:,i),lambda_1,param.optim.alpha,-1);
        tmp=lambdamize(param.optim.norm_bar,lambda,tmp);
        tmp=tmp.*(abs(tmp)>1e-16);
        subhess=subhess+tmp;
        tmp=hess_subroutine_smooth(x,eta,one_pic(:,i),lambda_2,param.optim.alpha,1);
        tmp=lambdamize(param.optim.norm_bar,lambda,tmp);       
        tmp=tmp.*(abs(tmp)>1e-16);
        subhess=subhess+tmp;
    end
end
              
       
            