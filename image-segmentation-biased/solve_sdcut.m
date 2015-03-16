function [cut, x_opt, x_real] = solve_sdcut(A, options)

    disp('-----------------------------------------------');
    disp('sdcut_solve_fast start ...');

    % scale A
    norm_A = sqrt(sum(sum(A.^2)));
    A = A ./ norm_A;

    % get B, b, m
    [B, B_non_sparse_part, b, c, a, gs, kappas, l_bbox, u_bbox, m, B_non_sparse_part_fast] ...
        = sdcut_pack_cons_data(A, options);

    % get start point
    u_init = calc_u_init(A, options, m);
       
    % solve dual
    tic
    [u_opt, iters] = solve_dual_lbfgsb_v3(u_init, A, B, b, l_bbox, u_bbox, ...
        options.sigma, B_non_sparse_part, options.lbfgsb_factr, options.lbfgsb_pgtol, options.lbfgsb_m, B_non_sparse_part_fast);
    t = toc;
    fprintf('lbfgsb take %.3f sec, %d iters (including line-search)\n', t, iters);

    % get optimal lifted primal X = x * x' 
    C_minus_opt = calc_c_minus(u_opt, A, B); 
    X_opt = (-0.5 / options.sigma) * C_minus_opt;

    % recover x from X = x * x'
    [x_opt, x_real] = sdcut_rand_hyper_rounding(X_opt, A, 10000, options, gs, kappas);    

    n = size(A, 1);
    cut = zeros(n, 2);
    cut(:,1) = x_opt == 1;
    cut(:,2) = x_opt == -1;

    disp('sdcut_solve_fast end');
    disp('-----------------------------------------------');
end

function u_init = calc_u_init(A, options, m)

    u_init = zeros(m,1);
    
end


function [u_opt, iters] = solve_dual_lbfgsb_v3(u_init, A, B, b, l_bbox, u_bbox, ...
    sigma, B_non_sparse_part, lbfgsb_factr, lbfgsb_pgtol, lbfgsb_m, B_non_sparse_part_fast)

    fcn = @(u) calc_dual_obj_grad_lbfgsb(u, A, B, b, sigma, B_non_sparse_part, B_non_sparse_part_fast);

    opts = struct('x0', u_init, 'maxIts', 1000, ...
        'factr', lbfgsb_factr, 'pgtol', lbfgsb_pgtol, 'm', lbfgsb_m);

    m = length(u_init);

    [u_opt, ~, info] = lbfgsb(fcn, l_bbox, u_bbox, opts);

    iters = info.totalIterations;

end


function [obj, grad] = calc_dual_obj_grad_lbfgsb(u, A, B, b, sigma, B_non_sparse_part, B_non_sparse_part_fast)

    %% calc C_minus
    persistent pre_V;
    persistent pre_D;

    num_neg_eigen = length(find(diag(pre_D) < 0));
    % 
    if isempty(pre_D) || num_neg_eigen > 100
        [C_minus, pre_V, pre_D] = calc_c_minus(u, A, B);
    else
        C = [];    

        eigen_opts.eigen_solver = 'eigs';
        eigen_opts.pre_V = pre_V;
        eigen_opts.pre_D = pre_D;
        eigen_opts.num_neg_eigen = num_neg_eigen + 2;
        eigen_opts.A_fun = @(x) calc_Cx(u, A, B, x, B_non_sparse_part, B_non_sparse_part_fast);

        [~, C_minus, pre_V, pre_D] = calc_pos_neg_part(C, eigen_opts); 
    end

    %% calc objective
    obj = calc_dual_obj(u, A, B, b, sigma, C_minus);

    %% calc gradient
    m = length(u);
    grad = zeros(m, 1);
    for ii = 1 : m
        grad(ii) = C_minus(:)' * B{ii}(:);
    end
    grad = (0.5 / sigma) * grad + b;

end


function Cx = calc_Cx(u, A, B, x, B_non_sparse_part, B_non_sparse_part_fast)

    n = size(A,1);


    % sparse part
    Cx = (sparse(1:n, 1:n, u(1:n)*B{1}(1,1), n, n) - A) * x;
    
    
    % structural part
    V  = B_non_sparse_part_fast.V;
    dd = B_non_sparse_part_fast.dd;
    Cx = Cx + V * (u(n+1:end) .* dd .* (V'*x));

end


function [C_minus, V, D] = calc_c_minus(u, A, B)

    m = length(B);

    C = -1 * A;
    for ii = 1 : m
        C = C + u(ii) * B{ii};
    end

    [~, C_minus, V, D] = calc_pos_neg_part(C, struct('eigen_solver', 'eig'));

end

function obj = calc_primal_obj(X, A, sigma)

    obj = -1 * A(:)' * X(:) + sigma * X(:)' * X(:);

end

function obj = calc_dual_obj(u, A, B, b, sigma, C_minus)

    if nargin == 5
        C_minus = calc_c_minus(u, A, B);
    end
        
    obj = (-0.25 / sigma) * sum(C_minus(:) .^ 2)  - u' * b;
    obj = -1 * obj;

end


function [X_plus, X_minus, V, D] = calc_pos_neg_part(X, options)


if nargin == 1
    options.eigen_solver = 'eig';
end


% make sure X is symmetric
% X = (X + X') / 2;

% eigenvectors and eigenvalues
switch options.eigen_solver
    case 'eig'
        [V, D] = eig(X);
    case 'eigs'
        max_d = -inf;
        k = options.num_neg_eigen;
        [n, K] = size(options.pre_V);
        eigs_opts.issym  = 1;
        eigs_opts.isreal = 1;  
        
        eigs_opts.v0 = options.pre_V * ones(K, 1);
        % eigs_opts.v0 = options.pre_V * rand(K, 1);
        % eigs_opts.v0 = options.pre_V(:,1);
            
        while max_d < 0
            eigs_opts.p = max(k, min(n, 2 * k + 10));
            fprintf(1, 'eigs start, k = %d ...\n', k);
            [V, D] = eigs(options.A_fun, n, k, 'sa', eigs_opts);
%             [V, D] = eigs(X, k, 'sa', eigs_opts);
            fprintf(1, 'eigs end\n');            
            
            max_d = max(diag(D));
            if max_d < 0
                min_d = min(diag(D));
                fprintf('max_d == %f < 0, min_d = %f, k = %d, n = %d\n', max_d, min_d, k, n);
                k = min(n, k*2);
            end
        end

    otherwise
        error('unknown eigen_solver: %s\n', eigen_solver);
end

dd = real( diag(D) );
idxs_minus = find(dd < 0);
n_minus = length(idxs_minus);
D_minus = sparse(1:n_minus, 1:n_minus, dd(idxs_minus), n_minus, n_minus, n_minus); 
V_minus = V(:, idxs_minus);
X_minus = (V_minus * D_minus) * V_minus';

fprintf('rank of negative part: %d, pos part: %d\n', n_minus, length(find(dd>0)));

if ~isempty(X)
    X_plus = X - X_minus;

    X_diff = max(max(abs(X - (X_plus + X_minus)))) / max(abs(X(:)));
    if max(X_diff(:)) > 0.001
        fprintf(1, 'X ~= X_plus + X_minus, diff: %f\n', max(X_diff(:)));
    end
else
    X_plus = [];
end

end

