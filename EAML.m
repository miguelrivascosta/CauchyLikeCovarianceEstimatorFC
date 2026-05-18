classdef EAML
    properties
        n
        nrf
        nrf2
        Gv
        Bm_dic
        Q
        M
        Av_grid
        theta_grid
        epsilon
    end
    methods
        function o = EAML(n,nrf)
            o.n = n;
            o.nrf = nrf;
            o.nrf2 = nrf*nrf;
            [o.Gv,o.Bm_dic,o.M] = o.gen_codebook();
            % o.theta_grid = deg2rad(theta0_deg:theta_step_deg:theta1_deg);
            o.theta_grid = linspace(-pi/2,pi/2,1500);
            o.Q = length(o.theta_grid);
            o.epsilon = 0.001;
            
            % Optimized grid pre-computation using pagemtimes
            SV = o.a_ula(pi*sin(o.theta_grid)); % (n, Q)
            o.Av_grid = pagemtimes(o.Gv, SV); % (nrf, Q, M)
        end

        function [Gv,Bm_dic,M] = gen_codebook(o)
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
            Gv = zeros(o.nrf,o.n,M);
	        for i = 1 : M
		        Bm_dic(:,:,i) = F(:,IDX(i,:)+1);
                Gv(:,:,i) = Bm_dic(:,:,i)';
	        end
        end

        function y = a_ula(o,psi)
            y = exp(1j*(0:o.n-1).'.*psi);
        end

        function hat_theta = estimate_doas(o,X,L)
            [~,K] = size(X);
            Km = K/o.M;
            % 1. Covariance Estimation
            Xm = reshape(X,[o.n,Km,o.M]);
            Ym = pagemtimes(o.Gv,Xm);
            Rv_hat = pagemtimes(Ym,pagectranspose(Ym))/Km;
            
            % 2. Whitening Transformation (Tv)
            inv_Rv_hat = pageinv(Rv_hat);
            k_scale = nthroot(Km,4);
            Tv = zeros(o.nrf,o.nrf,o.M);
            for v = 1 : o.M
                % Using sqrtm on pages. For small nrf, this is fast.
                Tv(:,:,v) = k_scale * sqrtm(inv_Rv_hat(:,:,v));
            end
            
            % Transformed manifolds and reference vector
            rhov_vec = sqrt(Km)*reshape(inv_Rv_hat,[o.M*o.nrf2,1]); 
            Av_tilde = pagemtimes(Tv,o.Av_grid);
            
            % 3. Vectorized Projection Component Pre-computation
            % Replaces the nested kron loops. A_cal_mat is (M*nrf2, Q)
            A_cal_mat = zeros(o.nrf2 * o.M, o.Q);
            for v = 1:o.M
                Av_v = squeeze(Av_tilde(:,:,v)); % (nrf, Q)
                % Vectorized equivalent of vec(a*a') for all columns
                block_v = reshape(Av_v, [o.nrf, 1, o.Q]) .* conj(reshape(Av_v, [1, o.nrf, o.Q]));
                A_cal_mat((v-1)*o.nrf2 + (1:o.nrf2), :) = reshape(block_v, [o.nrf2, o.Q]);
            end
            
            % Target vector y (vec_I) and its squared norm
            y = repmat(sqrt(Km)*reshape(eye(o.nrf),[o.nrf2,1]),[o.M,1]);
            yTy = y' * y;
            
            % 4. Multi-Signal Search with Successive Refinement
            idx_signals = zeros(1, L);
            
            % Initial Greedy Search (First Signal)
            loss = o.compute_loss_fast(A_cal_mat, rhov_vec, y, yTy, []);
            [~, idx_signals(1)] = min(loss);
            
            if L > 1
                % Successive search for signals 2 to L
                for l = 2:L
                    loss = o.compute_loss_fast(A_cal_mat, rhov_vec, y, yTy, idx_signals(1:l-1));
                    [~, idx_signals(l)] = min(loss);
                end
                
                % Iterative Refinement (Alternating Projection style)
                iter = 1;
                max_iter = 100;
                theta_hat_old = o.theta_grid(sort(idx_signals,'ascend'));
                while true
                    for l = 1:L
                        fixed_idx = idx_signals(1:L ~= l);
                        loss = o.compute_loss_fast(A_cal_mat, rhov_vec, y, yTy, fixed_idx);
                        [~, idx_signals(l)] = min(loss);
                    end
                    theta_hat = o.theta_grid(sort(idx_signals,'ascend'));
                    if all(abs(theta_hat_old - theta_hat) <= o.epsilon)
                        break;
                    end
                    theta_hat_old = theta_hat;
                    if max_iter == iter
                        break;
                    else
                        iter = iter + 1;
                    end
                end
            end
            hat_theta = o.theta_grid(sort(idx_signals,'ascend'));

        end

        function loss = compute_loss_fast(o, A_cal_mat, r, y, yTy, fixed_idx)
            % Optimized projection loss using the Schur Complement
            % Base correlations for fixed signals and rhov_vec (r)
            Q = size(A_cal_mat, 2);
            loss = zeros(1, Q);
            
            Af = A_cal_mat(:, fixed_idx);
            A_base = [Af, r];
            M_base = real(A_base' * A_base); % Real-part trick for complex LS
            D_base = real(A_base' * y);
            inv_M_base = inv(M_base);
            gain_base = D_base' * (inv_M_base * D_base);
            
            % Cross-correlations between the grid and fixed/reference components
            G_ab = real(A_cal_mat' * A_base); % (Q, length(fixed)+1)
            G_aa = sum(real(A_cal_mat .* conj(A_cal_mat)), 1).'; % (Q, 1)
            D_a = real(A_cal_mat' * y); % (Q, 1)
            
            mask = true(1, Q);
            mask(fixed_idx) = false;
            
            % Loop over the grid points to update the projection gain
            % This uses the Schur complement to solve (L+1)x(L+1) systems in scalar time
            D_base_transformed = inv_M_base * D_base;
            for q = find(mask)
                g_q = G_ab(q, :).';
                d_q = D_a(q);
                temp_v = inv_M_base * g_q;
                
                % Schur complement of the current candidate signal
                schur = G_aa(q) - g_q' * temp_v;
                % Update for the projection gain
                residual_d = d_q - g_q' * D_base_transformed;
                zq = residual_d / schur;
                
                % Loss = Total Energy - Projection Gain
                loss(q) = yTy - (gain_base + zq * residual_d);
            end
            loss(~mask) = NaN;
        end
    end
end