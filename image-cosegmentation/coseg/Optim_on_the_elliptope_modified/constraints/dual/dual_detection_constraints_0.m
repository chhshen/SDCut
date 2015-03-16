function subdual=dual_detection_constraints_0(x,param,eta)


global one_pic_bar
global one_pic


global thresh
global thresh2

global lambda


 global C
 
 subdual=zeros(size(C));
                                          
 for i=1:param.pic.pics
     if ~param.optim.smoothing
        [T2,T1]=dual_detection_constraints_subroutine(x,one_pic_bar(:,i),one_pic(:,i),thresh(i),thresh2(i),param);
        subdual=subdual+dual_detection_constraints_subroutine_2(one_pic_bar(:,i),one_pic(:,i),T1,T2); 
     else
        [T2,T1]=dual_detection_constraints_subroutine_smooth(x,one_pic_bar(:,i),one_pic(:,i),thresh(i),thresh2(i),param);
        subdual=subdual+dual_detection_constraints_subroutine_smooth_2(one_pic_bar(:,i),one_pic(:,i),T1,T2); 
     end    
 end
 
subdual=lambdamize(param.optim.norm,lambda,subdual);