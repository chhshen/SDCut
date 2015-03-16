function space=dual_detection_constraints_subroutine_smooth_2(one_pic_bar,one_pic,T1,T2)



I=find(one_pic);
J=find(one_pic==0);
I=repmat(I,1,size(J,1));
J=repmat(J',size(I,1),1);

space=sparse(I(:),J(:),1,size(one_pic,1),size(one_pic,1));

clear I J


space=sum(T2+T1)*(space);
space=space+space';

space=space+(diag((T2+T1).*one_pic)); 
