classdef PhaseIndependent
    properties
        n               % Number of antennas
        nrf             % Number of RF chains
        M               % Codebook size
    end
    
    methods
        function o = PhaseIndependent(n, nrf)
            % Constructor to initialize array and search parameters
            o.n = n;
            o.nrf = nrf;
            o.M = 2*o.n;
        end
        
        function hat_theta = estimate_doas(o, x, L)
            % Main method to estimate DOAs
            % x: Received signal matrix (N x T)
            % L: Number of sources
            
            % 1. Recover the Covariance Matrix Ry_hat using the shared logic
            Ry_hat = o.recover_covariance(x);
            hat_theta = rootmusic(Ry_hat,L);
        end
        
        function Ry_hat = recover_covariance(o, Y_dig)
            % Implementation of the Shmonin et al. recovery logic
            [N, size_data] = size(Y_dig);

            
            % Setup codebook and beams
            num_rep = size_data / (o.M / o.nrf);
            k = -o.n : o.n-1;
            psi_k = k * pi / o.n;
            B = 1/sqrt(N) * o.a_ula(psi_k);
            
            codebook_size = o.M / o.nrf;
            idx_mat = zeros(o.nrf, codebook_size);
            codebook = zeros(N, o.nrf, codebook_size);
            
            for i = 1 : codebook_size
                idx = i : 2*N/o.nrf : 2*N; 
                idx_mat(:,i) = idx.';
                codebook(:,:,i) = B(:,idx);
            end
            
            % Average power measurements
            avg_abs_y = zeros(2*N, 1);
            t = 1;
            for i = 1 : codebook_size
                Bm = codebook(:,:,i);
                for j = 1 : num_rep
                    idx = idx_mat(:, i);
                    % Power measurement for the beamspace components
                    avg_abs_y(idx) = avg_abs_y(idx) + abs(Bm' * Y_dig(:, t)).^2;
                    t = t + 1;
                end
            end
            avg_abs_y = avg_abs_y / num_rep;
            f_hat = reshape(avg_abs_y, [o.M, 1]);
            
            % Fourier series reconstruction of the covariance sequence
            q_vec = 1-N : N-1;
            exp_fourier_series = exp(-1j * pi * k.' * q_vec / o.n);
            rq = zeros(length(q_vec), 1);
            
            for ii = 1:length(q_vec)
                Dq = 1/(2*o.n) * sum(f_hat .* exp_fourier_series(:, ii), 1);
                rq(ii) = Dq / (N - abs(q_vec(ii)));
            end
            
            % Construct Toeplitz Covariance Matrix
            c = rq(N:-1:1);
            r = rq(N:1:end);
            Ry_hat = toeplitz(c, r);
        end
        
        function y = a_ula(o, psi)
            % Uniform Linear Array steering vector generator
            % psi can be a scalar or a row vector of spatial frequencies
            y = exp(1j * (0:o.n-1).' * psi);
        end
    end
end