classdef BSA
    properties
        N
        Nrf
        M
        codebook
        idx_mat;
        pseudo_inverse
    end
    
    methods
        function obj = BSA(N, Nrf)
            obj.N = N;
            obj.Nrf = Nrf;
            
            n_vec = (0:N-1);
            F = dftmtx(N)/sqrt(N);
            D = diag(exp(1j*pi/N*n_vec));
            F_hat = D*F;
            F_ext = zeros(N,2*N);
            F_ext(:,1:2:end) = F;
            F_ext(:,2:2:end) = F_hat;
            
            [~,idx] = sort(angle(F_ext(2,:)),'ascend');
            B = F_ext(:,idx);
            
            codebook_size = 2*N/Nrf;
            codebook = zeros(N,Nrf,codebook_size);
            idx_mat = zeros(Nrf,codebook_size);
            
            idx = 1:Nrf;
            for i = 1 : codebook_size
                idx = i : 2*N/Nrf : 2*N;
                idx_mat(:,i) = idx.';
                codebook(:,:,i) = B(:,idx);
                % idx = idx_mat(:,i) + Nrf;
            end
            
            C = [];
            for i = 1 : 2*N
                an = B(:,i);
                bn = kron(an,conj(an));
                C = [C;bn.'];
            end
            obj.pseudo_inverse = (C'*C+0.1*eye(power(N,2)))\C';
            obj.M = codebook_size;
            obj.idx_mat = idx_mat;
            obj.codebook = codebook;
            
        end

        function [R_hat,avg_abs_y] = reconstruction(obj,x,num_rep)
            avg_abs_y = zeros(2*obj.N,1);
            t = 1;
            for i = 1 : obj.M
                Bm = obj.codebook(:,:,i);
                for j = 1 : num_rep
                    idx = obj.idx_mat(:,i);
                    avg_abs_y(idx) = avg_abs_y(idx) + power(abs(Bm'*x(:,t)),2);
                    % avg_abs_y(idx) = avg_abs_y(idx) + diag(Bm'*R*Bm);
                    t = t + 1;
				end
				avg_abs_y(idx) = avg_abs_y(idx)/num_rep;
            end
            r_hat = obj.pseudo_inverse*avg_abs_y;
            R_hat = reshape(r_hat,[obj.N,obj.N]);
        end

        function hat_theta = estimate_doas(o,x,L)
            [N_dim,K] = size(x);

            % Ensure the input Y matches the expected operator dimensions
            if  N_dim ~= o.N
                error('Input matrix Y dimensions do not match initialized parameters.');
            end
            
            Km = K/o.M;
            % Ensure the ratio K/M is integer
            if  mod(Km,1)~=0
                error('The ratio K/M is not integer');
            end
            R_hat = o.reconstruction(x,Km);
            hat_theta = rootmusic(R_hat,L);
        end
    end
end