function subdual=dual_detection_constraints_mod3(x,npics,one_pic,...
                                              thresh,thresh2,lambda)


 global C
 
 subdual=zeros(size(C));
 
 n=size(x,1);
 
 for i=1:npics
            [space,T2,T1]=dual_detection_constraints_subroutine(x,one_pic(:,i),thresh(i),thresh2(i));
            
            T1=T1.*lambda(1+n+2*(i-1)*n:2*i*n);
            T2=T2.*lambda(1+2*(i-1)*n:n+2*(i-1)*n);
            
            %%%%%%%%%%%%%%%
            subdual=subdual+dual_detection_constraints_subroutine_2(space,one_pic(:,i),T1,T2);  
 end