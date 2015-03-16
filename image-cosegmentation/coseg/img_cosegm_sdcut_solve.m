%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Optimize the following sdp problem:  
% max   pobj(x) = <A, x*x'> + sigma*norm(x*x')^2                     
% s.t.  <B{i}, x*x'> = (or <=) b(i), i = 1, ..., m
%       x*x' is p.s.d.
% 
% - The corresponding dual problem is
% min   dobj(u) = 1/(4*sigma).*norm(C_minus)^2 + u'*b
% s.t.  u(i) >= 0, if i-th primal cons is '<=' 
% where C(u) = sum_i(u(i) .* B{i}) - A, C_minus(u) is negative part of C;
%       dgrd(u)(i) = 1/(2*sigma).*<C_minus(u), B{i}> + b(i)
% 
% - Relationship between primal and dual optimal vars:
% X_opt = -1/(2*sigma)*C_minus(u_opt)
%
% - Essential elements in struct 'data':
% n: num of primal vars
% m: num of primal cons
% A: mat(n,n)
% B: cell(1,m)
% b: vec(m,1)
% l_bbox: vec(m,1), l_bbox(i) = -inf, if i-th cons is '='; 
%                             = 0,    if i-th cons is '<='
% u_bbox: vec(m,1), u_bbox(i) = +inf, for any i 
% B_non_sparse_part: cell(m_non_sparse,1), see calc_Cx()
% fun_calc_Cx:   function handle, C(u)*x = func_calc_Cx(x, u, data)
% fun_calc_dobj: function handle, dobj(u) = func_calc_dobj(u, data, sigma, C_minus)
% fun_calc_dgrd: function handle, dgrd(u) = func_calc_dgrd(u, data, sigma, C_minus)
%
% - Format of large matrix (e.g. A, B, X, C, C_minus)
% 'full':            standard matlab 'full' matrix
% 'sparse':          standard matlab 'sparse' matrix
% 'part':            X_full = X.V*X.D*X.V'; X.V is full, X.D is sparse diagnal;
%                    cols of X.V are eigen-vectors, diag(D) are eigen-values. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cut, x_opt, scores, X_opt, u_opt, C_minus_opt, obj_val1, obj_val2] = img_cosegm_sdcut_solve(data, options)

    disp('-----------------------------------------------');
    disp('sdpcut_solve_fast start ...');
    
    
    % get B, b, l_bbox, u_bbox, m, etc ...
    %cons_data = feval(options.fun_pack_cons_data, options);
    
    % get start point
    if isempty(data.fun_calc_uinit)
        u_init = calc_u_init(data.A, data.m, options);
    else
        u_init = data.fun_calc_uinit(data, options);
    end
        
    % solve dual
    tic
    [u_opt, C_minus_opt, iters] = solve_dual_lbfgsb(u_init, data, options);
    t = toc;
    fprintf('solve_dual_lbfgsb takes %.3f sec, %d iters (including line-search)\n', t, iters);

    
    
    % get optimal lifted primal X = x * x'     
    X_opt = calc_X_opt(C_minus_opt, data, options);
  
    n = size(X_opt.V,1);
    diag_X     = sum(X_opt.V*X_opt.D.*X_opt.V, 2);
    inv_diag_X = 1 ./ sqrt(diag_X);
    X_opt.V = full(spdiags(inv_diag_X, 0, n,n) * X_opt.V);
        
    
    % recover x from X = x * x'
    [x_opt, scores] = rounding_ncut(X_opt, data, options);    

    [obj_val1, obj_val2] = check_correctness(X_opt, x_opt, u_opt, C_minus_opt, data, options);

    n = data.n;
    cut = zeros(n, 2);
    cut(:,1) = x_opt == 1;
    cut(:,2) = x_opt == -1;

    disp('sdpcut_solve_fast end');
    disp('-----------------------------------------------');
end

function X_opt = calc_X_opt(C_minus_opt, data, options)

switch data.outfmt
    case 'full'
        assert(strcmp(C_minus_opt.fmt, 'full'));
        X_opt.fmt = 'full';
        X_opt.mat = (-0.5 / options.sigma) * C_minus_opt.mat;
    case 'part'
        assert(strcmp(C_minus_opt.fmt, 'part'));
        X_opt.fmt = 'part';
        X_opt.D   = (-0.5 / options.sigma) * C_minus_opt.D;    
        X_opt.V   = C_minus_opt.V;
    otherwise
        error('invalid outfmt: %s\n', data.outfmt);
end

end


function u_init = calc_u_init(A, m, options)
    
    method = 'ncut';

    switch method
        case 'ncut'
            [~, x_init] = ncut_solve(A);
            n = size(A,1);
            u_init = [(A - 2 *(options.sigma)*eye(n)) * x_init ./ x_init; zeros(m-n,1)]; 
        case 'all-zero'
            u_init = zeros(m,1);
        case 'all-one'
            u_init = ones(m,1);
        case 'random'
            u_init = rand(m,1);
        otherwise 
            error('unknown method : %s\n', method);
    end
    
  
end

% function [u_opt, iters] = solve_dual_lbfgsb_v3(u_init, A, B, b, l_bbox, u_bbox, ...
%     sigma, B_non_sparse_part, lbfgsb_factr, lbfgsb_pgtol, lbfgsb_m)
function [u_opt, C_minus_opt, iters] = solve_dual_lbfgsb(u_init, data, options)
    global pre_V;
    global pre_D;
    global num_neg_eigen;
    global num_eig_iters;

    pre_V = [];
    pre_D = [];
    num_neg_eigen = data.num_neg_eigen_init;
    num_eig_iters = 0;

    fcn = @(u) calc_dual_obj_grad_lbfgsb(u, data, options);

    lbfgsb_opts = struct('x0',     u_init, ...
                         'maxIts', 10000, ...
                         'factr',  options.lbfgsb_factr, ... 
                         'pgtol',  options.lbfgsb_pgtol, ...
                         'm',      options.lbfgsb_m);

    [u_opt, ~, info] = lbfgsb(fcn, data.l_bbox, data.u_bbox, lbfgsb_opts);

    iters = info.totalIterations;
    
    C_minus_opt = calc_c_minus_new(u_opt, data, 'eigs', pre_V, pre_D, num_neg_eigen + 2);
    
end

function [C_minus, V, D] = calc_c_minus_new(u, data, method, pre_V, pre_D, num_neg_eigen)

    switch method
        case 'eig'
            C = calc_c(u, data.A, data.B); 
            %
            eigen_opts.eigen_solver = 'eig';
            eigen_opts.outfmt       = data.outfmt;
            %
            [~, C_minus, V, D] = calc_pos_neg_part(C, eigen_opts);
        case 'eigs'
            C = []; 
            %
            eigen_opts.eigen_solver = 'eigs';
            eigen_opts.n = data.n;
            eigen_opts.pre_V = pre_V;
            eigen_opts.pre_D = pre_D;
            eigen_opts.num_neg_eigen = num_neg_eigen;
            eigen_opts.outfmt = data.outfmt;
            if isempty(data.fun_calc_Cx)
                eigen_opts.A_fun = @(x) calc_Cx(u, data.A, data.B, x, data.B_non_sparse_part);
            else    
                eigen_opts.A_fun = @(x) data.fun_calc_Cx(C, x, u, data);
            end
            %
            [~, C_minus, V, D] = calc_pos_neg_part(C, eigen_opts); 
        otherwise
            error('unknown method: %s\n', method);
    end
end

%%
function [obj, grad] = calc_dual_obj_grad_lbfgsb(u, data, options)
    global pre_V;
    global pre_D;
    global num_neg_eigen;
    global num_eig_iters;

    A                 = data.A;
    B                 = data.B;
    b                 = data.b;
    sigma             = options.sigma;

    if 0
        [C_minus, pre_V, pre_D] = calc_c_minus_new(u, data, 'eig');
    else
        if isempty(num_neg_eigen) || num_neg_eigen > data.n * 0.75
            [C_minus, pre_V, pre_D] = calc_c_minus_new(u, data, 'eig');
        else
            if num_eig_iters > 40
                kk = 2;
            else
                kk = 5;
            end
            
            [C_minus, pre_V, pre_D] = calc_c_minus_new(u, data, 'eigs', ...
                                          pre_V, pre_D, num_neg_eigen + kk);
        end
    end
    num_neg_eigen = length(find(diag(pre_D) < 0));
    
    % calc objective
    if isempty(data.fun_calc_dobj)
        obj  = calc_dual_obj(u, A, B, b, sigma, C_minus.mat);
    else
        obj  = data.fun_calc_dobj(u, data, sigma, C_minus);
    end
    % calc gradient
    if isempty(data.fun_calc_dgrd)
        grad = calc_dual_grad(u, B, b, sigma, C_minus.mat);
    else
        grad = data.fun_calc_dgrd(u, data, sigma, C_minus);
    end
    
    num_eig_iters = num_eig_iters + 1;
    
end

% rule1: first n cons should be: X(i,i) == 1, for any i == 1 : n;
% rule2: sparse cons are arranged before structural cons
% rule3: structural cons can be expressed as: inner_prod(X, sign*(vec*vec')) <= b
function Cx = calc_Cx(u, A, B, x, B_non_sparse_part)

    n = size(A,1);
    m = length(B);
    m_non_sparse = length(B_non_sparse_part);

    % sparse part
    C_sparse = sparse(1:n, 1:n, u(1:n) * B{1}(1,1), n, n, n) - A;
    for ii = n+1 : m-m_non_sparse
        C_sparse = C_sparse + u(ii) * B{ii};
    end
    Cx_sparse = C_sparse * x;

    
    % structural part
    Cx_non_sparse = zeros(n,1);
    for ii = 1 : m_non_sparse
        vec    = B_non_sparse_part{ii}.vec;
        sign = B_non_sparse_part{ii}.sign;
        Cx_non_sparse = Cx_non_sparse + u(m-m_non_sparse+ii) * sign * (vec' * x) * vec; 
    end

    
    Cx = Cx_sparse + Cx_non_sparse;
end

function [C_minus, V, D] = calc_c_minus(u, A, B)

    C = calc_c(u, A, B); 
    
    eigen_opts.eigen_solver = 'eig';
    eigen_opts.outfmt       = 'full';
    
    [~, C_minus, V, D] = calc_pos_neg_part(C, eigen_opts);
    C_minus = C_minus.mat;

end

function C = calc_c(u, A, B)

    m = length(B);
    C = -1 * A;
    for ii = 1 : m
        C = C + u(ii) * B{ii};
    end
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

function grad = calc_dual_grad(u, B, b, sigma, C_minus)
    m = length(u);
    grad = zeros(m, 1);
    for ii = 1 : m
        grad(ii) = C_minus(:)' * B{ii}(:);
    end
    grad = (0.5 / sigma) * grad + b;
end

%%
function [x_opt, scores] = rounding_ncut(X, data, options)

    n = data.n;
    V = X.V;
    D = X.D;
    
    
%     diag_X     = sum(V*D.*V, 2);
%     inv_diag_X = 1 ./ sqrt(diag_X);
%     V = full(spdiags(inv_diag_X, 0, n,n) * V);

    % Assume W = 0.5 * (X + 1_nn)
    % Then Dg_vec = W*1_n = 0.5*(V*D*V'*1n + 1_nn*1n) = 0.5*V*D*sum(V',2) + 0.5*n
    %      Dg = diag(Dg_vec) --> Degree matrix 
    %      Dg_ha_vec = Dg_vec.^(-0.5)
    %      Dg_ha = diag(Dg_ha_vec) --> (Degree matrix).^(-0.5)
    %      Lsym = I - Dg_ha*W*Dg_ha = I - 0.5*Dg_ha*(X+1_nn)*Dg_ha
    %           = I - 0.5*Dg_ha*X*Dg_ha - 0.5*Dg_ha_vec*Dg_ha_vec'
    %      Lsym*x = x - 0.5*Dg_ha*V*D*V'*Dg_ha*x - 0.5*Dg_ha_vec*Dg_ha_vec'*x
    Dg_vec = 0.5*V*D*sum(V,1)' + 0.5*n;
    Dg_ha_vec = Dg_vec.^(-0.5);
    Dg_ha = diag(sparse(Dg_ha_vec));
    fun_calc_Lsym_x = @(x) x - 0.5*(Dg_ha*(V*(D*(V'*(Dg_ha*x))))) ...
                             - 0.5*(Dg_ha_vec*(Dg_ha_vec'*x));

    % ncut
    eigs_opts.issym  = 1;
    eigs_opts.isreal = 1;  
    fprintf(1, 'eigs start, k = 2 ...\n')
    [V,E] = eigs(fun_calc_Lsym_x, n, 2, 'sa', eigs_opts);
    fprintf(1, 'eigs end\n');
    z = V(:,2);


    x_opt = sign(z);
    x_opt(x_opt==0) = 1;
    scores = z;
end

%%
function [obj_val1, obj_val2] = check_correctness(X_opt, x_opt, u_opt, C_minus_opt, data, options)

    disp('----------------------------');
    disp('check_correctness start ...');

    %%
    if strcmp(X_opt.fmt, 'part')
        % check X_opt
        %data.fun_check_cons_X(X_opt, data, true);
        % check x_opt
        %data.fun_check_cons_x(x_opt, data, true);
        %check dual gap
        %pobj = data.fun_calc_pobj(X_opt, data, options.sigma);
        %dobj = data.fun_calc_dobj(u_opt, data, options.sigma, C_minus_opt) * -1;
        %fprintf(1, '- pobj == %f, dobj = %f, dual gap = %f\n', pobj, dobj, abs(pobj-dobj));
        % check bound's tightness
        obj_val1 = data.fun_calc_inner_AX(X_opt, data);
        obj_val2 = x_opt' * data.fun_calc_Ax(x_opt, data);
        fprintf(1, '- obj_val1 = %f\n', obj_val1);
        fprintf(1, '- obj_val2 = %f\n', obj_val2);
    end
end

function [X_plus, X_minus, V, D] = calc_pos_neg_part(X, options)
    % options.eigen_solver
    % options.pre_V;
    % options.pre_D;
    % options.A_fun;
    % options.k;

    % fprintf('options.eigen_solver: %s\n', options.eigen_solver);

    if nargin == 1
        options.eigen_solver = 'eig';
    end

    % eigenvectors and eigenvalues
    switch options.eigen_solver
        case 'eig'
            if strcmp(X.fmt, 'full')
                XX = (X.mat + X.mat') / 2;
                [V, D] = eig(XX);
            elseif strcmp(X.fmt, 'sparse')
                XX = full((X.S + X.S') / 2);
                [V, D] = eig(XX);
            else
                error('invalid X.fmt: %s\n', X.fmt);
            end

        case 'eigs'
            max_d = -inf;
            k = options.num_neg_eigen;
            n = options.n;

            eigs_opts.issym  = 1;
            eigs_opts.isreal = 1; 
            % eigs_opts.disp  = 2;
            % eigs_opts.maxit = 300;

            if ~isempty(options.pre_V)
                K = size(options.pre_V, 2);
                eigs_opts.v0 = options.pre_V * ones(K, 1);
                % eigs_opts.v0 = options.pre_V * rand(K, 1);
                % eigs_opts.v0 = options.pre_V(:,1);
            end

            while max_d < 0
                eigs_opts.p = max(k+1, min(n, 2 * k + 10));

                fprintf(1, 'eigs start, k = %d ...\n', k);
                %tic
                [V, D] = eigs(options.A_fun, n, k, 'sa', eigs_opts);
                % [V, D] = eigs(X, options.k, 'sa', eigs_opts);
                %t = toc;
                %fprintf(1, 'eigs end, %.3f sec eclipse\n', t);
                fprintf(1, 'eigs end\n');            

                max_d = max(diag(D));

                if max_d < 0
                    min_d = min(diag(D));
                    fprintf('max_d == %f < 0, min_d == %f, k = %d, n = %d\n', max_d, min_d, k, n);

                    if k >= n/2
                         tst_val = abs(max_d / min_d);
                         if tst_val > 0.001
                             error('abs(max_d / min_d) == %f > 0.001, k = %d, n = %d\n', tst_val, k, n);
                         else
                             break;
                         end
                    end                

                    k = min(n/2, k*2);
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

    fprintf('rank of negative part: %d\n', n_minus);

    % output X_minus
    if strcmp(options.outfmt, 'full')
        X_minus.fmt = 'full';
        X_minus.mat = (V_minus * D_minus) * V_minus';
    elseif strcmp(options.outfmt, 'part')
        X_minus.fmt = 'part';
        X_minus.V = V_minus;
        X_minus.D = D_minus;
    else
        error('invalid outfmt: %s\n', options.outfmt);
    end

    % output X_plus
    if ~isempty(X) && strcmp(options.outfmt, 'full')
        X_plus.fmt = 'full';
        if strcmp(X.fmt, 'full')
            X_plus.mat = X.mat - X_minus.mat;
            X_diff = max(max(abs(X.mat - (X_plus.mat + X_minus.mat)))) / max(abs(X.mat(:)));
        elseif strcmp(X.fmt, 'sparse')
            X_plus.mat = X.S - X_minus.mat;
            X_diff = max(max(abs(X.S - (X_plus.mat + X_minus.mat)))) / max(abs(X.S(:)));
        else
            error('invalid X.fmt: %s\n', X.fmt);
        end


        if max(X_diff(:)) > 0.001
            fprintf(1, 'X ~= X_plus + X_minus, diff: %f\n', max(X_diff(:)));
        end
    else
        X_plus = [];
    end
end





