function generate_kernel_laplacian(param)
    % COMPUTE KERNEL
    fprintf('COMPUTE KERNEL...')
    tic
    [xtilde,lW,param.optim.tab_lambda0,weights]     = main_kernel(param);
    t   = toc;
    fprintf('\n Elapsed time for the computation of the kernel is %f seconds\nDONE\n\n',t)

    [~, n_px] = size(xtilde);

    % Laplacian Matrix
    fprintf('COMPUTE LAPLACIAN MATRIX L...\n')
    tic
    L   = LaplacianMatrix_pixelwise(param,lW,n_px);
    t   = toc;
    fprintf(' Elapsed time for the computation of L is %f seconds\nDONE\n\n',t)

    % Chunks into SUPERPIXELS
    % deltas: indicator vectors n * n_pic, deltas(i,j) = 1 if ith fea point is in image j
    [deltas] = find_indices_image(lW,ones(sum(lW),1));

    
    dst_dir = [param.path.im,'kernel_laplacian/'];
    dst_file_path = [dst_dir, param.objet.type, '_kernel_laplacian.mat'];
    mkdir(dst_dir);
    save(dst_file_path, 'xtilde', 'lW', 'weights', 'n_px', 'L', 'deltas');
end