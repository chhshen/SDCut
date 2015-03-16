function [x_opt, x_real] = sdcut_rand_hyper_rounding(X, A, max_iter, options, gs, kappas)

[x_opt, x_real] = rand_hyper_rounding_srh(X, A, max_iter, options, gs, kappas);

end

%%
% SRH (Swept Random Hyperplanes) by K. Lang, NIPS2005. 
% 1. Randomly sample a vector u from N(0, X), and sort the nodes according to u;
% 2. Sweep a neighborhood of the cut sign(u), and find an optimal offset b such that 
%    the cut sign(u+b) has the largest obj value (x'*A*x).
% 3. Try several iters, filter out those trials violating constriants, 
%    and output the cut with largest obj value (x'*A*x).
% The following code is NOT efficient.
function [x_opt, x_real] = rand_hyper_rounding_srh(X, A, max_iter, options, gs, kappas)

    disp('----------------------------');
    fprintf(1, 'random hyperplane rounding start\n');

    same_pairs = options.same_pairs;
    diff_pairs = options.diff_pairs;
    
    n = size(X, 1);
    A1 = A - diag(diag(A));
    
    x_opt = zeros(n, 1);
    obj_val_opt = -inf;
    
    n_pass = 0;
    for ii = 1 : max_iter
        
        if mod(ii-1,10000) == 0
            us = mvnrnd(zeros(1,n), X, 10000)';
        end

        u = us(:,mod(ii-1,10000)+1);
       
        [~, order_flip_sign] = sort(u, 'descend');
        
        start_idx = max(1,   length(find(u>0)));
        end_idx   = min(n-1, length(find(u>0)));     

        x_tst = -1 * ones(n,1);
        x_tst(order_flip_sign(1:start_idx-1)) = 1;        
        obj_val_tst = x_tst' * A * x_tst;   
        
        x_opt_ii = -1 * ones(n,1);
        obj_val_opt_ii = -inf;
        
        for jj = start_idx : end_idx
                        
            idx_flip_sign = order_flip_sign(jj);
            x_tst(idx_flip_sign) = +1;
            
            obj_val_tst = obj_val_tst + A1(idx_flip_sign, :) * x_tst * 4;
            
            check_flag = check_same_diff_cons(x_tst, same_pairs, diff_pairs);
            if check_flag == 0
                continue;
            end

            check_flag = check_bias_cons(x_tst, gs, kappas);
            if check_flag == 0
                continue;
            end
            
            if(obj_val_tst > obj_val_opt_ii)
                obj_val_opt_ii = obj_val_tst;
                x_opt_ii = x_tst;                
            end
            
        end
        
        if ~isempty(find(x_opt_ii > 0, 1))
            n_pass = n_pass + 1;
        end
        
        if obj_val_opt_ii > obj_val_opt
            obj_val_opt = obj_val_opt_ii;
            x_opt = x_opt_ii;
            x_real = u;
        end

    end
    
    disp('random hyperplane rounding end');
    disp('----------------------------');
end

%%
function check_flag = check_same_diff_cons(x, same_pairs, diff_pairs)
        
    check_flag = 1;
    for idx = 1 : size(same_pairs, 2)
        i = same_pairs(1, idx);
        j = same_pairs(2, idx);
        if x(i) ~= x(j)
            check_flag = 0;
        end
    end

    for idx = 1 : size(diff_pairs, 2)
        i = diff_pairs(1, idx);
        j = diff_pairs(2, idx);
        if x(i) == x(j)
            check_flag = 0;
        end
    end
    
end

function check_flag = check_bias_cons(x, gs, kappas)
    check_flag = 1;
    for idx = 1 : length(kappas)
        g = gs(:,idx);
        kappa = kappas(idx);
        if abs(g' * x) < kappa * 0.99
            check_flag = 0;
        end
    end
end


