function [B, B_non_sparse_part, b, c, a, gs, kappas, l_bbox, u_bbox, m, B_non_sparse_part_fast] ...
    = sdcut_pack_cons_data(A, options)

n = size(A, 1);

same_pairs = options.same_pairs;
diff_pairs = options.diff_pairs;

bias_group_pos  = options.bias_group_pos;
bias_group_neg  = options.bias_group_neg;
bias_kappa_same = options.bias_kappa_same;
bias_kappa_diff = options.bias_kappa_diff;


%% count #constraints
% #cons (sparse) for diag(X) = e
m_cons_diag = n; 

% #cons (sparse) for same/diff 
m_cons_same = size(same_pairs, 2); 
m_cons_diff = size(diff_pairs, 2); 

% #cons (structural) for bias group
m_cons_bias = 0;
if ~isempty(bias_group_pos)
    m_cons_bias = m_cons_bias + 1;
end
if ~isempty(bias_group_neg)
    m_cons_bias = m_cons_bias + 1;
end
if ~isempty(bias_group_pos) && ~isempty(bias_group_neg)
    m_cons_bias = m_cons_bias + 1;
end

% #cons (structural) for balancing
switch options.balance_cons_type  
    case 'equpartition',     m_cons_bala = 1; 
    case 'non-equpartition', m_cons_bala = 2;
    case 'normalized',       m_cons_bala = 2;
    otherwise, error('unknown options.balance_cons_type: %s\n', options.balance_cons_type);
end

% number of primal constraints
m = m_cons_diag + m_cons_same + m_cons_diff + m_cons_bias + m_cons_bala; 
m_non_sparse = m_cons_bias + m_cons_bala;


%% construct B and b
B = cell(m,1);
B_non_sparse_part = cell(m_non_sparse, 1);
b = zeros(m, 1);
l_bbox = zeros(m, 1);
u_bbox = zeros(m, 1);

% cons (sparse) for diag(X) = e
for ii = 1 : m_cons_diag
    B{ii} = sparse(ii,ii,1,n,n);
    b(ii) = 1; 
    l_bbox(ii) = -inf;
    u_bbox(ii) = inf;
end

% cons (sparse) for same/diff 
scale_cons_same = 1 / sqrt(2);
scale_cons_diff = scale_cons_same;
for ii = 1 : m_cons_same
    i = same_pairs(1, ii);
    j = same_pairs(2, ii);
    B{ii+m_cons_diag} = scale_cons_same * sparse([i,j],[j,i],[1,1],n,n);
    b(ii+m_cons_diag) = scale_cons_same * 2;
    l_bbox(ii+m_cons_diag) = -inf;
    u_bbox(ii+m_cons_diag) = inf;
end
for ii = 1 : m_cons_diff
    i = diff_pairs(1, ii);
    j = diff_pairs(2, ii);
    B{ii+m_cons_diag+m_cons_same} = scale_cons_diff * sparse([i,j],[j,i],[1,1],n,n);
    b(ii+m_cons_diag+m_cons_same) = scale_cons_diff * -2;
    l_bbox(ii+m_cons_diag+m_cons_same) = -inf;
    u_bbox(ii+m_cons_diag+m_cons_same) = inf;
end

% cons (structural) for bias group
gs     = zeros(n, m_cons_bias);
kappas = zeros(1, m_cons_bias);
B_idx = m_cons_diag+m_cons_same+m_cons_diff;
idx = 0;

dd = sum(A,2);
D_inv = diag(dd.^(-1));
P = D_inv * A;

if ~isempty(bias_group_pos)
    B_idx = B_idx + 1;
    idx = idx + 1;
    T = bias_group_pos;
    kappa = bias_kappa_same;

    g = zeros(n,1);
    g(T) = 1;
    g = g'*P;
    g = g';

    scale_ratio = 1 / norm(g);
    g = scale_ratio * g;
    kappa = sum(g) * kappa;

    B{B_idx} = -1 * (g * g');
    b(B_idx) = -1 * (kappa^2);
    l_bbox(B_idx) = 0;
    u_bbox(B_idx) = inf;
    B_non_sparse_part{idx}.vec = g;
    B_non_sparse_part{idx}.sign = -1;
    gs(:,idx) = g;
    kappas(idx) = kappa;
end

if ~isempty(bias_group_neg)
    B_idx = B_idx + 1;
    idx = idx + 1;
    T = bias_group_neg;
    kappa = bias_kappa_same;

    g = zeros(n,1);
    g(T) = 1;
    g = g'*P;
    g = g';

    scale_ratio = 1 / norm(g);
    g = scale_ratio * g;
    kappa = sum(g) * kappa;

    B{B_idx} = -1 * (g * g');
    b(B_idx) = -1 * (kappa^2);
    l_bbox(B_idx) = 0;
    u_bbox(B_idx) = inf;
    B_non_sparse_part{idx}.vec = g;
    B_non_sparse_part{idx}.sign = -1;
    gs(:,idx) = g;
    kappas(idx) = kappa;
end

if ~isempty(bias_group_pos) && ~isempty(bias_group_neg)
    B_idx = B_idx + 1;
    idx = idx + 1;
    T_pos = bias_group_pos;
    T_neg = bias_group_neg;
    kappa = bias_kappa_diff;

%     A1 = A - eye(n);
    A1 = A - diag(diag(A));
    d = sum(A1, 2);
    Vol_T_pos = sum(d(T_pos));
    Vol_T_neg = sum(d(T_neg));

    g = zeros(n,1);
    g(T_pos) = +1 / length(T_pos);
    g(T_neg) = -1 / length(T_neg);
    g = g'*P;
    g = g';

    scale_ratio = 1 / norm(g);
    g = scale_ratio * g;
    kappa = sum(abs(g)) * kappa;

    B{B_idx} = -1 * (g * g');
    b(B_idx) = -1 * (kappa^2);
    l_bbox(B_idx) = 0;
    u_bbox(B_idx) = inf;
    B_non_sparse_part{idx}.vec = g;
    B_non_sparse_part{idx}.sign = -1;
    gs(:,idx) = g;
    kappas(idx) = kappa;   
end


% cons (structural) for balancing
switch options.balance_cons_type 

    % a2^2 <= (c'*x)^2 = < c*c', X > <= a1^2
    case 'non-equpartition'
        scale_cons_balance = 1 / sqrt(n);
        
        c    = scale_cons_balance * ones(n,1);
        a.a1 = scale_cons_balance * options.a1 * n;
        a.a2 = scale_cons_balance * options.a2 * n;
        
        B{end-1} = c * c';
        b(end-1) = a.a1^2;
        l_bbox(end-1) = 0;
        u_bbox(end-1) = inf;
        B_non_sparse_part{end-1}.vec = c;
        B_non_sparse_part{end-1}.sign = 1;
        
        B{end} = -1 * (c * c');
        b(end) = -1 * (a.a2^2);
        l_bbox(end) = 0;
        u_bbox(end) = inf;
        B_non_sparse_part{end}.vec = c;
        B_non_sparse_part{end}.sign = -1;
    
    otherwise
        error('unknown options.balance_cons_type: %s\n', options.balance_cons_type);
end


m_non_sparse = length(B_non_sparse_part); 
B_non_sparse_part_fast.V = [];
B_non_sparse_part_fast.dd = [];

for ii = 1 : m_non_sparse
    vec = B_non_sparse_part{ii}.vec;
    sign = B_non_sparse_part{ii}.sign; 
    B_non_sparse_part_fast.V = [B_non_sparse_part_fast.V, vec];
    B_non_sparse_part_fast.dd = [B_non_sparse_part_fast.dd; sign];
end



end
