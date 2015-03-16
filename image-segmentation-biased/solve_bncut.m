%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An implementation of the biased normalized cut framework described in the paper:
% 
% Biased Normalized Cuts,
% Subhransu Maji, Nisheeth Vishnoi and Jitendra Malik,
% In Proceedings, CVPR 2011, Colorado Springs.
% http://www.cs.berkeley.edu/~smaji/projects/biasedNcuts/
%
% The original code was modified by Peng Wang

function [cut, x_opt, x_real] = solve_bncut(A, options)

nvec   = options.biased_ncut_nvec;
bias_group_pos = options.bias_group_pos;

[EigVal, EigVect, D] = affinity2evec(A,nvec);

% calc s_T
n = size(A, 1);
bw = zeros(n,1);
for i = 1 : length(bias_group_pos)
    bw(bias_group_pos(i)) = 1;
end
s = D*bw(:);

% calc gamma
gamma = -1 * mean(EigVal(2:end));

% calc linear combination
wts = zeros(length(EigVal),1);
for i = 2:length(EigVal) % ignore the all 1's vector
    wts(i) = (EigVect(:,i)'*s)/(EigVal(i) - gamma);
end
WEigVect = EigVect*wts;

x_real = WEigVect; 
x_opt  = sign(x_real);
x_opt(x_opt==0) = 1;

cut = zeros(n, 2);
cut(:,1) = x_opt == 1;
cut(:,2) = x_opt == -1;

norm_A = sqrt(sum(sum(A.^2)));
A = A ./ norm_A;
fprintf(1, 'x''*A*x: %f\n', x_opt'*A*x_opt);

end

function [EigVal, EigVect, D] = affinity2evec(A,nvec)

[wx, wy] = size(A);
x = 1 : wx;
S = full(sum(A, 1));
D = sparse(x, x, S, wx, wy);
clear S x;

opts.issym=1;
opts.isreal = 1;
% opts.disp=2;
[EigVect, EVal] = eigs(D - A, D, nvec, 'sm',opts);
clear W opts;

EigVal = diag(EVal);

EigVal(1:end) = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);

% txo=orig_sz(1); tyo=orig_sz(2); 
% vect = zeros(txo*tyo, nvec);
% for v = 1 : nvec,
%     rszEigVect = imresize(reshape(EigVect(:, v), [ty tx])',[txo,tyo]);
%     vect(:,v)  = rszEigVect(:);
% end

end



