function [A, norm_A] = compute_A(param, xtilde, L, lW_px, lW, weights, P)

if nargout == 1
    A = P'*P-sum(P)'*sum(P)/sum(lW_px)-(xtilde*P)' * xtilde*P+param.optim.lapWght.*P'*L*P;
else
    
    % A = S + V*D*V';       
    S = sparse(P'*P + param.optim.lapWght.*P'*L*P);
    V = full([(xtilde*P)', sum(P)' / sqrt(sum(lW_px))]); 
    D = sparse(-1 * eye(size(V,2)));

    A.fmt = 'sparse+struct';
    A.S = S;
    A.V = V; 
    A.D = D;
    

    % norm(S+V*D*V')^2 = Tr((S + V*D*V')(S + V*D*V')) 
    %                  = Tr(S*S + S*V*D*V' + V*D*V'*S + V*D*V'*V*D*V')
    %                  = norm(S)^2 + 2*sum((S*V*D).*V)(:)) + norm(V'*V*D)^2
    norm_A = sqrt( sum(sum(S.^2)) + 2*sum(sum((S*V*D).*V)) + sum(sum((V'*V*D).^2)) );

end

end

% function [A, Ax_diag_PtP, Ax_sum_P, Ax_sum_lWpx, Ax_xtildeP, Ax_lapWght, Ax_PtLP] ...
%          = compute_A(param, xtilde, L, lW_px, lW, weights, P)
% 
% % sqrtD*sqrtD' <=> P*1_n*1_n'*P'
% 
% % le eye(p) remplacer par diag(sqrtD) ?
% 
% %A   = eye(p) - sqrtD * sqrtD' /n - xtilde' * xtilde + lambdaAJ.*Laj;
% 
% 
% A = P'*P-sum(P)'*sum(P)/sum(lW_px)-(xtilde*P)' * xtilde*P+param.optim.lapWght.*P'*L*P;
% 
% 
% % C*x =   diag_PtP .* x 
% %       - sum_P' * ((sum_P * x) / sum_lWpx) 
% %       - xtildeP' * (xtildeP * x) 
% %       + lapWght * (PtLP * x);
% Ax_diag_PtP = diag(P'*P);
% Ax_sum_P    = sum(P);
% Ax_sum_lWpx = sum(lW_px);
% Ax_xtildeP  = xtilde*P;
% Ax_lapWght  = param.optim.lapWght;
% Ax_PtLP     = sparse(P'*L*P);
% 
% A_data.sparse = sparse(P'*P + param.optim.lapWght.*P'*L*P);
% A_data.struct.V = [(xtilde*P)', sum(P)' / sqrt(sum(lW_px))]; 
% A_data.struct.D = sparse(-1 * eye(size(A_data.struct.V,2)));
% 
% A_tst = A_data.sparse + A_data.struct.V * A_data.struct.D * A_data.struct.V';
% 
% keyboard
% 
% return
% 
% sum_lW_px_i=0;
% sum_lW_i=0;
% for i=1:param.pic.np_considered+param.pic.nd_considered
% 
%     subC=eye(lW_px(i))-ones(lW_px(i))/sum(lW_px);
% 
%     subC=subC-...
%         xtilde(: , sum_lW_px_i+1:sum_lW_px_i+lW_px(i))'...
%         *xtilde(: , sum_lW_px_i+1:sum_lW_px_i+lW_px(i));
% 
%     subC=subC+param.optim.lapWght.*L(sum_lW_px_i+1:sum_lW_px_i+lW_px(i),sum_lW_px_i+1:sum_lW_px_i+lW_px(i));
% 
% 
%     %rajouter les poids :
%     if i<=param.pic.np_considered
%         wght_i=1./sqrt(weights(1));
%     else
%         wght_i=1./sqrt(weights(2));
%     end
%     
%     subC=wght_i.*subC.*wght_i;
%     
%     C(sum_lW_i+1:sum_lW_i+lW(i),sum_lW_i+1:sum_lW_i+lW(i))=...
%         P(sum_lW_px_i+1:sum_lW_px_i+lW_px(i),sum_lW_i+1:sum_lW_i+lW(i))'...
%         *subC...
%         *P(sum_lW_px_i+1:sum_lW_px_i+lW_px(i),sum_lW_i+1:sum_lW_i+lW(i));
% 
%     clear subC
% 
%     sum_lW_px_j=sum_lW_px_i+lW_px(i);
%     sum_lW_j=sum_lW_i+lW(i);
% 
% 
%     for j=i+1:param.pic.np_considered+param.pic.nd_considered
%         subC=-ones(lW_px(i),lW_px(j))/sum(lW_px);
% 
%         subC=subC-...
%             xtilde( :, sum_lW_px_i+1:sum_lW_px_i+lW_px(i))'...
%             *xtilde( : , sum_lW_px_j+1:sum_lW_px_j+lW_px(j));
% 
%         subC=subC+param.optim.lapWght.*L(sum_lW_px_i+1:sum_lW_px_i+lW_px(i) , sum_lW_px_j+1:sum_lW_px_j+lW_px(j));
% 
%         
%         
%         if j<=param.pic.np_considered
%             wght_j=1./sqrt(weights(1));
%         else
%             wght_j=1./sqrt(weights(2));
%         end
% 
%         subC=wght_j*subC*wght_i;
% 
%         try
%             C(sum_lW_i+1:sum_lW_i+lW(i),sum_lW_j+1:sum_lW_j+lW(j))=...
%                 P(sum_lW_px_i+1:sum_lW_px_i+lW_px(i),sum_lW_i+1:sum_lW_i+lW(i))'...
%                 *subC...
%                 *P(sum_lW_px_j+1:sum_lW_px_j+lW_px(j),sum_lW_j+1:sum_lW_j+lW(j));
%         catch
%             keyboard
%         end
%         clear subC
% 
%         C(sum_lW_j+1:sum_lW_j+lW(j),sum_lW_i+1:sum_lW_i+lW(i))=C(sum_lW_i+1:sum_lW_i+lW(i),sum_lW_j+1:sum_lW_j+lW(j))';
% 
%         sum_lW_px_j=sum_lW_px_j+lW_px(j);
%         sum_lW_j=sum_lW_j+lW(j);
% 
%     end
%     sum_lW_px_i=sum_lW_px_i+lW_px(i);
%     sum_lW_i=sum_lW_i+lW(i);
% end
% 
% 
% % 
% % 
% % A   =diag(sqrtD.*sqrtD); 
% % 
% % 
% % 
% % 
% % A = A - sqrtD * sqrtD' /n ;
% % 
% % A = A - xtilde' * xtilde;
% % 
% % A= A + lambdaAJ.*Laj;