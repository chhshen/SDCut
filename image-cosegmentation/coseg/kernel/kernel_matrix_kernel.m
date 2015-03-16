function xtilde = kernel_matrix_kernel(x,df,n)

        K = x;

        if ~isempty(df);
            % df is given
            mmax = round(min(1000,4*df));   % maximum size for the incomplete Cholesky
            [G,P,m] = icd_full(K,n*1e-8,min(n,mmax)); % perform largest ICD to get lambda
            I = P(1:m);
            [temp,Pi] = sort(P);
            G = G(Pi,1:m);

            xtilde = G';
            xtilde = xtilde - repmat(mean(xtilde,2),1,size(xtilde,2));
            % [u,e] = eig(xtilde * xtilde');
            [u,e,v] = svd(xtilde,'econ');
            e=e.^2;
            ind =find( real(diag(e))/n > 1e-10 );
            if df >= length(ind),
                % df is too big
                lambda = 1e-10;
            else
                lambda = df_to_lambda_eig(df,real(diag(e(ind,ind)))/n);
                lambda = 1/lambda;
            end
            ind = find( real(diag(e)) > n*lambda*1e-2);
            q = u(:,ind);
            xtilde = q' * xtilde;
            C = xtilde * xtilde' + n * lambda * eye(size(xtilde,1));
            R = chol(C);
            R = inv(R)';
            xtilde = R * xtilde;

            
        else
            % lambda is given
            [G,P,m] = icd_full(K,n*lambda*1e-2,n);
            I = P(1:m);
            [temp,Pi] = sort(P);
            G = G(Pi,1:m);
            G = G';
            G = G - repmat(mean(G,2),1,size(G,2));
            C = G * G' + n * lambda * eye(size(G,1));
            try
                R = chol(C);
            catch
                % augment lambda if not positive enough and error is produced
                R = chol(C + norm(G,'fro')^2 * 1e-10 * eye(size(G,1)));
                warning('lambda was augmented');
            end
            R = inv(R)';
            xtilde = R * G;
            

            
        end