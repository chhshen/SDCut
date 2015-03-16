function [R]=constraints_subsubroutine(x,V,lambda,signe,param)



id_im=V*ones(1,size(x,2));%I_{im}
sum_x=sum(x.*id_im);% \bar{Y} dans le rapport
T=(x*sum_x');%Y*\bar{Y}'

if param.optim.smoothing==0
    R=x_plus(signe.*(T-lambda),param.optim.epsilon);
else
    R=signe.*(T-lambda);
end
