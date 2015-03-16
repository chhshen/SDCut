function lambda = df_to_lambda_eig(df,e);
n = length(e);
if df>n+1e-12, lambda = NaN; return; end
if df>n-1e-12, lambda = Inf; return; end


% % % binary search

a = 0;
b = .5 * ( df + n );
b = 1/min(e) * b/n / ( 1 - b/n);

i=1;
while (b-a)/(b+a)>1e-12
    i=i+1;
    c=(a+b)/2;
    if sum(  c * e./ ( 1 + c*e ) ) > df, b = c; else a=c; end
if (i>1000), keyboard; end
end
lambda =c;

% lambda = 1/mean(e);
% optparam.kmax    = 20;
% optparam.display = 1;
% optparam.tol     = 1e-20;
% optparam.tol_df  = 1e-20;
% [lambda,exit_parameters] = minimize_newton(lambda,@cost_df_to_lambda,optparam,e,df);
%     