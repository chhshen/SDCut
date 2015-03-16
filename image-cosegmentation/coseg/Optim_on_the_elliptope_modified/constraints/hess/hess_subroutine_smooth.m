function subhess=hess_subroutine_smooth(x,eta,one_pic,lambda_1,alpha,signe)


sum_x=sum_x_over_V(x,one_pic);

sum_eta=sum_x_over_V(eta,one_pic);

R=(x*sum_x');

R=signe*(R-lambda_1);

R=alpha*R;


[Rp,Rn]=seperate_signe(R);


%[V]=subroutine_logExp(Rp,Rn);

[U]=function_sigmExp(Rp,Rn);

%si au carre :
%T=V.*U;
% sinon :
T=U;

%subhess=signe*2*alpha*(one_pic*(T'*eta)+(T)*sum_eta);
subhess=signe*alpha*(one_pic*(T'*eta)+(T)*sum_eta);

%si carre :
%R=signe*(U.^2.+U.*(exp(-Rp)./(exp(-Rp)+exp(Rn)).*V);

R=signe*U.*(exp(-Rp)./(exp(-Rp)+exp(Rn)));

R_x=x_over_V(x,R);
R_eta=x_over_V(eta,R);

% subhess=subhess...
%         +2*alpha*alpha*(...
%         one_pic*( (x*sum_eta'+eta*sum_x')'*R_x)...
%         +R_eta*(sum_x'*sum_x)...
%         +R_x*(sum_eta'*sum_x)...
%         );
subhess=subhess...
        +alpha*alpha*(...
        one_pic*( (x*sum_eta'+eta*sum_x')'*R_x)...
        +R_eta*(sum_x'*sum_x)...
        +R_x*(sum_eta'*sum_x)...
        );
    
    
subhess=subhess./(alpha);





