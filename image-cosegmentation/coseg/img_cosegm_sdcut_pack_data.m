function data = img_cosegm_sdcut_pack_data(param)

%%% notation: 
% im=>images, 
% px=>subsampled-pixels, also feature-points
% sp=>superpixel
%
%%% vars: 
% bad_indice: 1     x n_sp0,  indice(max:n_sp) of bad sps 
% P:          n_px  x n_sp1,  P(i,j) = 1, if i-th px in j-th good-sp; 
%                                    = 0, otherwise
% deltas:     n_sp1 x n_im,   deltas(i,j) = #pxs in i-th good-sp, if this sp in j-th im
%                                         = 0, otherwise
% A:          n_sp1 x n_sp1,  affinity matrix used in objective function of SDP
% xtilde:     d     x n_px,   internal used for quick computation
% lW:         n_im  x 1,      #good-sps in each im     
% lW_px:      n_im  x 1,      #pxs in each im
% weights:
%
% n_im:        #images (ims)
% n_px:        #subsampled-pixels = #feature-points (pxs)
% n_sp:        #superpixels (sps), n_sp = n_sp0 + n_sp1
% n_sp0:       #bad--sps, which do not contain any px
% n_sp1:       #good-sps, which contain at least one px
% lW_sp:       n_im  x 1,      #sps in each im
% lW_sp0:      n_im  x 1,      #bad-sps in each im
% good_indice: 1 x n_sp1, indice(max:n_sp) of good sps 
% sp_maps:     cell(1, n_im), superpixel maps for each frame


%% pack function handles
if 1
data.outfmt = 'part'; % the formats of X, C, C_minus are 'part'. 
data.fun_calc_uinit     = @img_cosegm_calc_uinit;
data.fun_calc_Cx        = @img_cosegm_calc_Cx;
data.fun_calc_dobj      = @img_cosegm_calc_dobj;
data.fun_calc_dgrd      = @img_cosegm_calc_dgrd;
data.fun_calc_pobj      = @img_cosegm_calc_pobj;
data.fun_calc_inner_AX  = @img_cosegm_calc_inner_AX;
data.fun_calc_Ax        = @img_cosegm_calc_Ax;
data.fun_check_cons_x   = @img_cosegm_check_cons_x;
data.fun_check_cons_X   = @img_cosegm_check_cons_X;
data.num_neg_eigen_init = 10;
end


%% pack A data
kernel_laplacian_file = [param.path.im, 'kernel_laplacian/', ...
                         param.objet.type, '_kernel_laplacian.mat'];
load(kernel_laplacian_file);

% Compute the Discriminative clustering matrix A
fprintf('COMPUTE DISCRIMINATIVE CLUSTERING MATRIX A...\n')
tic
lW_px               = lW;
[P,lW,bad_indice,sp_maps]   = chunk_norris(param);
% [A, Ax_diag_PtP, Ax_sum_P, Ax_sum_lWpx, Ax_xtildeP, Ax_lapWght, Ax_PtLP] ...
%                             = compute_A(param,xtilde,L,lW_px,lW,weights,P);
[A, norm_A] = compute_A(param,xtilde,L,lW_px,lW,weights,P);
                        
deltas              = P'*deltas;



% trA = norm(A);%trace(A);
% A   = A/trA;
% % in our case, A is to maximized, instead of minimized.
% A   = -1 * A; 
A.S = (-1 /norm_A) * A.S;
A.D = (-1 /norm_A) * A.D;

% pack aux_data
data.A           = A;
data.deltas      = sparse(deltas);
data.n_im        = length(lW);
data.n_px        = n_px;
data.n_sp0       = length(bad_indice);
data.n_sp1       = size(P,2);
data.n_sp        = data.n_sp0 + data.n_sp1;
data.good_indice = setdiff(1:data.n_sp, bad_indice);
data.bad_indice  = bad_indice;
data.sp_maps     = sp_maps;
data.P           = P;
% data.xtilde      = xtilde;
data.lW          = lW;
data.lW_px       = lW_px;
% data.Ax_diag_PtP = Ax_diag_PtP; 
% data.Ax_sum_P    = Ax_sum_P;
% data.Ax_sum_lWpx = Ax_sum_lWpx; 
% data.Ax_xtildeP  = Ax_xtildeP; 
% data.Ax_lapWght  = Ax_lapWght; 
% data.Ax_PtLP     = Ax_PtLP;
% data.Ax_trA      = trA;



data.lW_sp       = zeros(data.n_im,1);
data.lW_sp0      = zeros(data.n_im,1);
for ii = 1 : data.n_im
    data.lW_sp(ii)  = length(unique(sp_maps{ii}));
    data.lW_sp0(ii) = data.lW_sp(ii) - lW(ii);
end

data.indice_sp_im = cell(1, data.n_im);
idx_start = 1;
for ii = 1 : data.n_im
    idx_end = idx_start + data.lW_sp(ii) - 1;
    data.indice_sp_im{ii} = idx_start : idx_end;
    idx_start = idx_end + 1;
end


%% pack cons data
n    = data.n_sp1;
p    = size(deltas, 2);
bala_cons_type = param.bala_cons_type;

assert(n == size(deltas,1));

% count #constraints
m_cons_diag = n;
switch bala_cons_type
    case 1
        m_cons_bala = 2 * n * p;
    case 2
        m_cons_bala = p;
    otherwise
        error('invalid bala_cons_type: %d\n', bala_cons_type);
end
m = m_cons_diag + m_cons_bala;
m_non_sparse = 0;


% construct B and b
B = cell(m,1);
B_non_sparse_part = cell(m_non_sparse, 1);
b = zeros(m, 1);
l_bbox = zeros(m, 1);
u_bbox = zeros(m, 1);

% cons (sparse) for diag(X) = e
coni = 1;
for ii = 1 : m_cons_diag
    clear B_i;
    B_i.fmt = 'sparse-one-elem';
    B_i.S = [ii,ii,1];%sparse(ii,ii,1,n,n,1);%%
    
    
    B{coni} = B_i;
    b(coni) = 1; 
    l_bbox(coni) = -inf;
    u_bbox(coni) = inf;
    coni = coni + 1;
    
end


% cons (sparse) for balance
if bala_cons_type == 1
    lamda0 = param.lamda0;
    lamda1 = param.lamda1;
    n_is = sum(deltas);
    deltas_scale = zeros(p,1);
    deltas_scaled = zeros(n,p);
    for ii = 1 : p
        delta_i = deltas(:,ii);
        deltas_scale(ii) = 1 / sqrt(sum(delta_i.^2));
        deltas_scaled(:,ii) = delta_i * deltas_scale(ii);
    end

    for ii = 1 : p
        n_i = n_is(ii);
        delta_i = deltas(:,ii);
        b1 = (2*lamda1-1)*n_i;
        %b0 = (1-2*lamda0)*n_i;
        for jj = 1 : n
            %B{coni} = sparse(n,n);
            %B{coni}(jj,:) = delta_i * deltas_scale(ii);
            B{coni} = [1, ii,jj];
            
            
            b(coni) = b1 * deltas_scale(ii);
            l_bbox(coni) = 0;
            u_bbox(coni) = inf;
            coni = coni + 1;
        end
    end

    for ii = 1 : p
        n_i = n_is(ii);
        delta_i = deltas(:,ii);
        %b1 = (2*lamda1-1)*n_i;
        b0 = (1-2*lamda0)*n_i;
        for jj = 1 : n
            %B{coni} = sparse(n,n);
            %B{coni}(jj,:) = -1 * delta_i * deltas_scale(ii);
            B{coni} = [-1, ii,jj];
            
            b(coni) = b0 * deltas_scale(ii);
            l_bbox(coni) = 0;
            u_bbox(coni) = inf;
            coni = coni + 1;
        end
    end
    data.deltas_scale = deltas_scale;
    data.lamda0 = lamda0;
    data.lamda1 = lamda1;
    data.deltas_scaled = deltas_scaled;
end

if bala_cons_type == 2
    lamda = param.lamda;
    n_is = sum(deltas);
    deltas_scale  = zeros(p,1);
    deltas_scaled = zeros(n,p);
    for ii = 1 : p
        n_i = n_is(ii);
        delta_i = deltas(:,ii);
        deltas_scale(ii) = 1 / sqrt(sum(delta_i.^2));%1 / (lamda * n_i);%% % %
        deltas_scaled(:,ii) = delta_i * deltas_scale(ii);
        
        clear B_i;
        B_i.fmt = 'struct';
        B_i.V = deltas_scale(ii) * delta_i;
        B_i.D = 1;
        
        % B{coni} = delta_i * delta_i' * (deltas_scale(ii)^2);
        B{coni} = B_i;       
        b(coni) = (lamda * n_i * deltas_scale(ii))^2;
        l_bbox(coni) = 0;
        u_bbox(coni) = inf;
        coni = coni + 1;
        
    end
    data.deltas_scale  = deltas_scale;
    data.deltas_scaled = sparse(deltas_scaled);
    data.lamda = lamda;
end


data.bala_cons_type = param.bala_cons_type;
data.B = B;
data.b = b;
data.B_non_sparse_part = B_non_sparse_part;
data.l_bbox = l_bbox;
data.u_bbox = u_bbox;
data.m = m;
data.n = n;



end

function dobj = img_cosegm_calc_dobj(u, data, sigma, C_minus)
    if ~strcmp(C_minus.fmt, 'part')
        error('invalid fmt: %s\n', C_minus.fmt);
    end
    
    dd = diag(C_minus.D);
    
    % norm(C_minus)^2 = Tr(C_minus*C_minus') = sum(dd.^2)
    dobj = (1/(4*sigma)) * sum(dd.^2) + u' * data.b;
end

function dgrd = img_cosegm_calc_dgrd(u, data, sigma, C_minus)
    if ~strcmp(C_minus.fmt, 'part')
        error('invalid fmt: %s\n', C_minus.fmt);
    end

    
    deltas_scaled = data.deltas_scaled;
    b  = data.b;
    m  = data.m;
    V  = C_minus.V;
    D  = C_minus.D;
    dd = diag(D);
    
    if data.bala_cons_type == 2
        %%% calc <C_minus, B{i}>, i=1:m
        % dgrd(1:n) = diag(V*D*V')
        % dgrd(n+i) = <V*D*V', del_i*del_i'> = (del_i'*V)*D*(V'*del_i), i = 1:p
        if ~isempty(dd)
            dgrd = [V.^2 * dd; (deltas_scaled'*V).^2 * dd];
        else
            dgrd = zeros(m,1);
        end
        dgrd = (1/(2*sigma)) * dgrd + b;
        
    elseif data.bala_cons_type == 1
        %%% calc <C_minus, B{i}>, i=1:m
        % dgrd(1:n) = diag(V*D*V')
        % dgrd(n+i*n+j)     = <V*D*V', [zeros(j-1,n);    del_i'; zeros(n-j,n)]> 
        %    =     (V*D*V'*dels)(j,i), i = 1:p, j = 1:n
        % dgrd(n+n*p+i*n+j) = <V*D*V', [zeros(j-1,n); -1*del_i'; zeros(n-j,n)]> 
        %    = -1* (V*D*V'*dels)(j,i), i = 1:p, j = 1:n
        if ~isempty(dd)
            dgrd1 = V.^2 * dd;
            dgrd2 = 2 * V*(D*(V'*deltas_scaled)); %%%%%%%%%%%%%%%%%%%%%%%
            dgrd3 = dgrd2 * -1;
            dgrd = [dgrd1; dgrd2(:); dgrd3(:)];
        else
            dgrd = zeros(m,1);
        end
        %%% 
        dgrd = (1/(2*sigma)) * dgrd + b;
    else
        error('invalid data.bala_cons_type: %s\n', data.bala_cons_type);
    end
end

function pobj= img_cosegm_calc_pobj(X, data, sigma)
    dd = diag(X.D);

    pobj = -1 * img_cosegm_calc_inner_AX(X, data);    
    pobj = pobj + sigma * sum(dd.^2);
end

function val = img_cosegm_calc_inner_AX(X, data)
    if ~strcmp(X.fmt, 'part')
        error('invalid fmt: %s\n', X.fmt);
    end   

    V = X.V;
    D = X.D;
    dd = diag(D);
    
    val = 0;
    for ii = 1 : length(dd)
        d = dd(ii);
        v = V(:,ii);
        val = val + d * v' * img_cosegm_calc_Ax(v, data);
    end
end


function Ax = img_cosegm_calc_Ax(x, data)

      S = data.A.S;
      V = data.A.V;
      D = data.A.D;

      Ax = S*x;
      Ax = Ax + (V*(D*(V'*x))); 
end


%%
function Cx = img_cosegm_calc_Cx(C, x, u, data)

% calc Ax
Ax = img_cosegm_calc_Ax(x, data);
 
% calc u(i)*B{i}*x
if data.bala_cons_type == 1
    p = data.n_im;
    dels = data.deltas_scaled;
    
    uBx_1 = u(1:n) .* x;
    
    u_mat2 = reshape(u(n+1    :n+n*p),   n, p);
    uBx_2 = u_mat2 * (dels' * x) + dels * (u_mat2' * x);
    
    u_mat3 = reshape(u(n+n*p+1:n+2*n*p), n, p);
    uBx_3 = -1 * ( u_mat3 * (dels' * x) + dels * (u_mat3' * x) );
    
    uBx = uBx_1 + uBx_2 + uBx_3;

elseif data.bala_cons_type == 2
    n = data.n;
    p = data.n_im;
    deltas_scaled = data.deltas_scaled;
    
    uBx_1 = u(1:n) .* x;
    uBx_2 = deltas_scaled * ((deltas_scaled' * x) .* u(n+1:n+p));   
    
    uBx = uBx_1 + uBx_2;
       
else
    error('invalid bala_cons_type: %d\n', data.bala_cons_type);
end

% 
Cx = uBx - Ax; 



end

%%
function u_init = img_cosegm_calc_uinit(data, options)

method = 'ncut';

A = data.A;
m = data.m;
n = data.n;


switch method
    case 'ncut'
        if 1
            n  = data.n;
            S  = A.S;
            V  = A.V;
            dd = diag(A.D);
            sigma = options.sigma;
            
            % Assume L = -A = Degree - W,
            % Then, Degree = diag(L) + diag(W) = diag(-A) + eye,
            % L_rw = Degree_inv * L = -1 * Degree_inv * A
            Degree_inv = diag( sparse( ( -1*(diag(S)+(V.*V)*dd) + 1 ) .^-1) );
            fun_calc_Lrw_x = @(x) -1*Degree_inv*img_cosegm_calc_Ax(x, data);
            
            % ncut
            eigs_opts.issym  = 1;
            eigs_opts.isreal = 1;  
            fprintf(1, 'eigs start, k = 2 ...\n');
            [V,E] = eigs(fun_calc_Lrw_x, n, 2, 'sa', eigs_opts);
            fprintf(1, 'eigs end\n');
            z = V(:,2);
            
            x_init = z;
            u_init = (img_cosegm_calc_Ax(x_init, data) - 2*sigma*x_init) ./ x_init;
            u_init = [u_init; zeros(m-n,1)];
        end
    case 'all-zero'
        u_init = zeros(m,1);
    case 'all-one'
        u_init = ones(m,1);
        norm_u = norm(u_init);
        u_init = u_init / norm_u;
        u_init(n+1:m) = 0;
    case 'random'
          u_init = rand(m,1);
    otherwise 
        error('unknown method : %s\n', method);
end

end
























