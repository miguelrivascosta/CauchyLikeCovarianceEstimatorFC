classdef CLGLS < handle
    properties
        n
		nrf
		Bm_dic
		M
		IDX
        Psi
        Sigmacl
        Omega
        Pf
        Phi

        alphal_dic
        alpha_mat
        
        P00
        P10
        ep
        Phi_sparse
        v
        
        Pf_sparse

        PhiPf

        ak
        bk

        p
        q
    end
    methods
		function o = CLGLS(n,nrf)
			o.n = n;
			o.nrf = nrf;
			[o.Bm_dic,o.IDX,o.M] = o.gen_codebook();
            [o.Psi,o.Sigmacl,o.Pf,o.Phi,o.Omega] = o.gen_linear_operators();
            o.alphal_dic = o.gen_alphal_dic(o.n,o.nrf);
            o.alpha_mat = o.gen_alpha_mat(o.n,o.nrf);
            
            o.p = 2*o.n+o.M-1;
            o.q = 2*(o.M*(o.nrf-1)-o.n)+1;
            o.Pf_sparse = sparse(o.Pf);
            o.P00 = sparse(o.Pf(1:o.p,1:end-1));
            o.P10 = sparse(o.Pf(o.p+1:end,1:end-1));
            o.ep = sparse([zeros(o.p-1,1);1]);
            v = zeros(2*o.n-1,1);
            v(2:2:end) = -1;
            o.v = sparse(v);
            o.Phi_sparse = sparse([eye(2*n-1);o.v.']);
            o.PhiPf = sparse(o.Phi.'*o.Pf.');
            [o.ak,o.bk] = o.gen_ak_bk(o.n);


        end
		
        function alpha_mat = gen_alpha_mat(~,n,nrf)
	        [V_mesh, U_mesh] = meshgrid(1:nrf, 1:nrf);
	        alpha_mat = (2j/n) ./ (1 - exp(1j*2*pi/n * (V_mesh - U_mesh)));
	        alpha_mat(logical(eye(nrf))) = 2;
        end

        function alphal_dic = gen_alphal_dic(~,n,nrf)
	        [V_mesh, U_mesh] = meshgrid(1:nrf, 1:nrf);
	        alphal_dic = -(2j/n) ./ (1 - exp(1j*2*pi/n * (V_mesh - U_mesh)));
	        alphal_dic(logical(eye(nrf))) = 0;
        end
		function y = a_ula(o,psi)
			y = exp(1j*(0:o.n-1).'.*psi);
        end

        function [R_hat,r_hat] = reconstruction(o,X,Km)
            Xm = reshape(X,[o.n,Km,o.M]);

            whitening_mat = zeros(o.nrf*o.nrf,o.nrf*o.nrf,o.M);
            inv_sm_vec = zeros(o.nrf*o.nrf,1,o.M);
            for m = 1 : o.M
                Bm = o.Bm_dic(:,:,m);
                X = Xm(:,:,m);
                Ym = Bm'*X;
                Sm = Ym*Ym'/Km;
                inv_Sm = eye(o.nrf)/Sm;
                inv_sm_vec(:,:,m) = inv_Sm(:);
                whitening_mat(:,:,m) = kron(inv_Sm.',inv_Sm);
            end
            aux0 = sum(pagemtimes(pagemtimes(pagectranspose(o.Lm_mat),whitening_mat),o.Lm_mat),3);
            inv_aux0 = eye(2*o.n-1)/aux0;
            aux1 = sum(pagemtimes(pagectranspose(o.Lm_mat),inv_sm_vec),3);
            r_hat = inv_aux0*aux1;
			aux = [real(r_hat(1));r_hat(2:2:end)+1j*r_hat(3:2:end)];
			R_hat = toeplitz(aux,conj(aux));

        end

        function [Bm_dic,IDX,M] = gen_codebook(o)
	        F = 1/sqrt(o.n)*exp(1j*2*pi/o.n*(0:o.n-1).'*(0:o.n-1));
	        
	        if o.nrf==o.n
		        M = 1;
	        else
		        M = ceil(o.n/(o.nrf-1));
	        end
	        
	        IDX = zeros(M,o.nrf);
	        IDX(1,:) = 0:o.nrf-1;
	        for i = 1 : M-1
		        IDX(i+1,:) = mod(IDX(i,:) + o.nrf-1,o.n);
	        end
	        
	        Bm_dic = zeros(o.n,o.nrf,M);
	        for i = 1 : M
		        Bm_dic(:,:,i) = F(:,IDX(i,:)+1); 
	        end
        end

        function [Psi,Sigmacl,Pf,Phi,Omega] = gen_linear_operators(o)
            Sigmacl = o.gen_Sigmacl(o.n,o.nrf);
            Phi = o.gen_Phi(o.n);
            Omega = o.gen_Omega(o.n);
            Pf = o.gen_Pf(o.n,o.nrf);
            Psi = kron(eye(o.M),Sigmacl)*Pf*Phi*Omega;
        end

        function Phi = gen_Phi(~,N)
	        v = zeros(2*N-1,1);
	        v(2:2:end) = -1;
	        Phi = [eye(2*N-1);v'];
        end
        
        function Omega = gen_Omega(~,N)
	        u_vec = 0:2*N-2;
	        q_vec = -N+1:N-1;
	        Omega = zeros(2*N-1);
            for u = 1 : length(u_vec)
                for q = 1 : length(q_vec)
                    Omega(u,q) = Omega_uq(u_vec(u),q_vec(q),N);
                end
            end

            function y = Omega_uq(u,q,N)
	            if mod(u,2) == 0
		            y = 1/2*(1-abs(q)/N)*exp(-1j*pi/N*u*q);
	            else
		            y = sin(pi*abs(q)/N)*exp(-1j*pi/N*u*q);
	            end
            end
        end
        
        function Sigma = gen_Sigmacl(~,n,nrf)
            alpha_uv = zeros(nrf,nrf);
            for u = 1 : nrf
                for v = 1 : nrf
                    if u~=v
                        alpha_uv(u,v) = 2j/n*1/(1-exp(1j*2*pi/n*(v-u)));
                    end
                end
            end
            Sigma = zeros(nrf*nrf,2*nrf-1);
            for l = 1 : nrf
                D = zeros(nrf);
                D(l,l) = 1;
                Sigma(:,2*l-1) = 2*D(:);
            end
            
            for l = 1 : nrf-1
                D = zeros(nrf);
                for u = 1 : l
                    for v = l+1 : nrf
                        Euv = zeros(nrf);
                        Euv(u,v) = 1;
        
                        D = D + alpha_uv(u,v)*Euv + conj(alpha_uv(u,v))*Euv.';
                    end
                end
                Sigma(:,2*l) = D(:);
            end
        end
        
        function Pf = gen_Pf(~,n, nrf)
            if n == nrf
                M_local = 1;
            else
                M_local = ceil(n / (nrf - 1));
            end
            rows_per_block = 2*nrf - 1;
            total_cols = 2*n;
            
            % 1. Pre-allocate the entire matrix for speed
            Pf = zeros(M_local * rows_per_block, total_cols);
            
            % 2. Define the base indices for the first block (Pf0)
            % Pf0 has 1s at (i, i) for i = 1 to 2*nrf-1
            base_row_indices = (1:rows_per_block)';
            base_col_indices = 1:rows_per_block;
            
            % 3. Fill the blocks using cyclic indexing
            shift_step = 2 * (nrf - 1);
            
            for m = 0 : M_local-1
                % Calculate the shift for this block
                current_shift = m * shift_step;
                
                % Calculate shifted column indices with wrap-around (1-based indexing)
                % This replaces Pf0 * PI^k
                shifted_cols = mod(base_col_indices + current_shift - 1, total_cols) + 1;
                
                % Calculate the vertical position in the large Pf matrix
                row_offset = m * rows_per_block;
                
                % Linear indexing to place 1s in the correct spots
                % We are placing a 1 at (row_offset + i, shifted_cols(i))
                target_rows = row_offset + base_row_indices;
                linear_indices = sub2ind(size(Pf), target_rows, shifted_cols');
                
                Pf(linear_indices) = 1;
            end
        end
        
        function hat_r = estimate_hatr(o,x,Km)
            x_reshaped = reshape(x,[o.n,Km,o.M]);
            inv_hat_Rep = zeros(o.M*o.nrf*o.nrf,o.M*o.nrf*o.nrf);
            hatp = zeros(o.nrf*o.nrf,o.M);

            ym = pagemtimes(pagectranspose(o.Bm_dic),x_reshaped);
            hat_Sm_dic = pagemtimes(ym,pagectranspose(ym))/Km;
            inv_hat_Sm_dic = pageinv(hat_Sm_dic);
            for m = 1 : o.M
                hat_Sm = hat_Sm_dic(:,:,m);
                inv_hat_Sm = inv_hat_Sm_dic(:,:,m);
                
                hatp(:,m) = hat_Sm(:);

                inv_hat_Repm = 1/Km*kron(inv_hat_Sm.',inv_hat_Sm);
                
                idx = (m-1)*o.nrf*o.nrf+1:m*o.nrf*o.nrf; 
                inv_hat_Rep(idx,idx) = inv_hat_Repm;
            end
            hatp = hatp(:);
            hat_r = (o.Psi'*inv_hat_Rep*o.Psi)\(o.Psi'*inv_hat_Rep*hatp);
            hat_r = [real(hat_r(o.n));hat_r(o.n+1:end)];
        end
        
        function hat_r = estimate_hatr_fast(o,x,Km)
            x_reshaped = reshape(x,[o.n,Km,o.M]);

            ym = pagemtimes(pagectranspose(o.Bm_dic),x_reshaped);
            hat_Sm_dic = pagemtimes(ym,pagectranspose(ym))/Km;
            hat_invSm_dic = pageinv(hat_Sm_dic);
            %% Stage I
            first_order_trace = zeros(o.M*(2*o.nrf-1),1);
            second_order_trace = zeros(o.M*(2*o.nrf-1));
            for m = 1 : o.M
               
                idx = (m-1)*(2*o.nrf-1)+1:m*(2*o.nrf-1);
                first_order_trace(idx,1) = fast_trace0(o.nrf,o.alphal_dic,hat_invSm_dic(:,:,m));
                second_order_trace(idx,idx) = fast_trace1(o.nrf,o.alpha_mat,o.alphal_dic,hat_invSm_dic(:,:,m));
            end
            
            %% Stage II
            M00 = sparse(second_order_trace(1:o.p,1:o.p));
            M01 = sparse(second_order_trace(1:o.p,o.p+1:end));
            M10 = sparse(second_order_trace(o.p+1:end,1:o.p));
            M11 = sparse(second_order_trace(o.p+1:end,o.p+1:end));
            Q = o.P00.'*M00*o.P00;
            U = [o.v o.P00.'*M00*o.ep+o.P10.'*M10*o.ep o.P10.' o.P00.'*M01];
            %C = [ [M00(p,p), 1; 1, 0], zeros(2,2*q); zeros(2*q,2), [M11, eye(q); eye(q), zeros(q,q)] ]; 
            C_inv = [ [0, 1; 1, -M00(o.p,o.p)], zeros(2, 2*o.q);zeros(2*o.q, 2), [zeros(o.q,o.q), eye(o.q); eye(o.q), -M11] ];

            %% Stage III
            hat_scl = eval_woodbury(Q,U,C_inv,o.PhiPf*first_order_trace);
            
            %%% Stage IV
            hat_r = inv_Omega(o.n,o.ak,o.bk,hat_scl);
            hat_r = [real(hat_r(1));hat_r(2:end)];
            function y = fast_trace0(nrf,alphal_dic,inv_Sm)
	            y = zeros(2*nrf,1);
	            y(1:2:end) = diag(inv_Sm);
	            y(2:2:end) = cumsum(sum(real(conj(alphal_dic) .* inv_Sm), 1).');
	            y = 2*y(1:end-1);
            end

            function y = fast_trace1(nrf,alpha_mat,alphal_dic,inv_Sm)
                y = zeros(2*nrf);
                l_idx = zeros(2*nrf-1,1);
                l_idx(1:2:end) = 1:nrf;
                l_idx(2:2:end) = 1:1:nrf-1;
            
                % Wl = zeros(nrf);
                Wl_old = zeros(nrf);
	            dot_prods = alphal_dic'*inv_Sm;
            
                for i = 1 : 2*nrf-1
                    l = l_idx(i);
                
                    if mod(i,2) == 1
                        Wl = 2*(inv_Sm(:,l)*inv_Sm(l,:));
                    elseif i == 2
    		            term1 = inv_Sm(:, 1) * (alpha_mat(1, :) * inv_Sm);
    		            term2 = (inv_Sm * alpha_mat(:, 1)) * inv_Sm(1, :);
    		            
    		            Wl_old = term1 - term2;
                        Wl = Wl_old;
                    else
                        aux = inv_Sm(:,l)*dot_prods(l,:);
                        Wl = Wl_old + aux + aux';
                        Wl_old = Wl;
                    end
                    
		            y(1:2:end,i) = diag(Wl);
		            y(2:2:end,i) = cumsum(sum(real(conj(alphal_dic) .* Wl), 1).');
                end
                y = 2*y(1:end-1,1:end-1);
            end
        
            function y_fast = eval_woodbury(Q,U,C_inv,x)
                a = Q \ x;
                B = Q \ U;
                K = C_inv + U' * B;
                y_fast = a - B * (K \ (U' * a));
            end
        
            function y = inv_Omega(n,ak,bk,scl)
	            X = ifft([scl;-sum(scl(2:2:end))]);
	            y = ak.*X(1:n)+bk.*X(n+1:end);
            end
        end

        function [ak,bk] = gen_ak_bk(~,n)
	        ak = zeros(n,1);
	        bk = zeros(n,1);
	        for k = 0 : n-1
		        ak(k+1) = (eval_we(n,k-n)+eval_wo(n,k-n))/(eval_we(n,k)*eval_wo(n,k-n)+eval_wo(n,k)*eval_we(n,k-n));
		        bk(k+1) = -(eval_we(n,k-n)-eval_wo(n,k-n))/(eval_we(n,k)*eval_wo(n,k-n)+eval_wo(n,k)*eval_we(n,k-n));
            end

            function we = eval_we(n,q)
	            we = 1/2*(1-abs(q)/n);
            end
            
            
            function wo = eval_wo(n,q)
	            wo = sin(pi/n*abs(q));
            end
        end

        function hat_theta = estimate_doas(o,x,L,options)
            arguments
                o
                x
                L
                options.mode (1,1) string {mustBeMember(options.mode, ["fast", "direct"])} = "fast"
            end
            % x: N x K  matrix
            % L: Number of sources to estimate
            
            [N_dim, K] = size(x);
            Km = K/o.M;
            % Ensure the input Y matches the expected operator dimensions
            if  N_dim ~= o.n
                error('Input matrix Y dimensions do not match initialized parameters.');
            end

            % Ensure the ratio K/M is integer
            if  mod(Km,1)~=0
                error('The ratio K/M is not integer');
            end
            
            if options.mode == "direct"
                hat_r = o.estimate_hatr(x,Km);
            else
                hat_r = o.estimate_hatr_fast(x,Km);
            end
            hat_R = toeplitz(hat_r,conj(hat_r));
    
            hat_theta = sort(rootmusic(hat_R,L),'ascend');
        end
    end
end
