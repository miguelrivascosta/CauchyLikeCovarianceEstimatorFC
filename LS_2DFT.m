classdef LS_2DFT
    properties
        N
        Nrf
        M
        idx_mat
        idx_diagonals
        codebook
        diag_constants
        F_base
        avg_mat
        complex_constant_mat
        pseudo_inverse
        max_l
        A
    end
    
    methods
        function obj = LS_2DFT(N, Nrf)
            obj.N = N;
            obj.Nrf = Nrf;
            obj.F_base = obj.dft_centered();
            
            if Nrf == N 
                M = 1;
            else
                M = ceil(N/(Nrf-1));
            end
            obj.M = M;

            idx_mat = zeros(Nrf,obj.M);
            codebook = zeros(N,Nrf,obj.M);
            idx_mat(:,1) = 0:Nrf-1;
            codebook(:,:,1) = obj.F_base(:,idx_mat(:,1)+1);
        
            for i = 1 : obj.M - 1
                idx_mat(:,i+1) = mod(idx_mat(:,i) + Nrf-1,N);
                codebook(:,:,i+1) = obj.F_base(:,idx_mat(:,i+1)+1);
                idx = idx_mat(:,i) + Nrf;
            end

            obj.idx_mat = idx_mat;
            obj.codebook = codebook;

            A = zeros(M*Nrf*(Nrf+1)/2,2*N-1);
            ctd = 1;
            for m = 0 : M - 1
                idx = idx_mat(:,m+1);
                for l = 0 : Nrf - 1
                    row = idx(1:end-l);
                    col = mod(row + l,N);
                    for z = 0 : length(row)-1
                        A(ctd,:) = obj.get_ecuation(row(z+1),col(z+1));
                        ctd = ctd + 1;
                    end
                end
            end
            
            obj.pseudo_inverse = (A'*A)\A';
            obj.A = A;
        end

        function y = row_diag(obj,m)
            y = zeros(1,2*obj.N-1);
            y(1) = 1;
            k_vec = 1:obj.N-1;
        
            y(2:obj.N) = (1 - k_vec/obj.N).*exp(1j*pi*k_vec).*exp(-1j*2*pi*m*k_vec/obj.N);
            y(obj.N+1:end) = conj(y(2:obj.N));
        end
        
        function y = row_off_diag(obj,m,n)
            y = zeros(1,2*obj.N-1);
            k_vec = 1:obj.N-1;
            y(2:obj.N) = exp(1j*pi*k_vec).*(exp(-2j*pi/obj.N*m*k_vec)-exp(-2j*pi/obj.N*n*k_vec));
            y(obj.N+1:end) = -conj(y(2:obj.N));
            y = 1/obj.N*1/(1-exp(1j*2*pi/obj.N*(n-m)))*y;
        end

        function y = get_ecuation(obj,m,n)
            if m == n
                y = obj.row_diag(m);
            else
                y = obj.row_off_diag(m,n);
            end
        end

        function R_hat = cauchy_like_recovery(obj, x, num_rep)            
            t = 1;
            vector = [];
            for i = 0 : obj.M - 1
                B = obj.codebook(:,:,i+1);
                Sm = zeros(obj.Nrf,obj.Nrf);
                for j = 0 : num_rep -1
                    y = B'*x(:,t);
                    Sm = Sm + y*y';
                    t = t + 1;
                end
                Sm = Sm / num_rep;
                for j = 0 : obj.Nrf-1
                    vector = [vector;diag(Sm,j)];
                end
            end
            
            
            r_hat = obj.pseudo_inverse*vector;
            
            r_hat(1) = real(r_hat(1));
            r_hat = r_hat(1:obj.N);
            R_hat = toeplitz(r_hat,conj(r_hat));
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
            R_hat = o.cauchy_like_recovery(x,Km);
            hat_theta = rootmusic(R_hat,L);
        end


        function R_hat = rec_mat(obj,zero_diagonal,first_diagonal)
            R_hat = zeros(obj.N,obj.N);
            for i = 0 : obj.N - 1
                for j = i : obj.N - 1
                    if (i == j)
                        R_hat(i+1,j+1) = zero_diagonal(j+1);
                    elseif (i+1 == j)
                        R_hat(i+1,j+1) = first_diagonal(j);
                        R_hat(j+1,i+1) = conj(first_diagonal(j));
                    else
                        idx = (i:j-1)+1;
                        R_hat(i+1,j+1) = obj.complex_constant_mat(i+1,j+1)*sum(first_diagonal(idx));
                        R_hat(j+1,i+1) = conj(R_hat(i+1,j+1));
                    end
                end
            end
        end
        
        function Fc = dft_centered(obj)
            n_vec = (0:obj.N-1);
            m_vec = (0:obj.N-1);
            Fc = 1/sqrt(obj.N)*exp(1j*2*pi/obj.N*(m_vec-obj.N/2).*n_vec.');
        end

    end
end