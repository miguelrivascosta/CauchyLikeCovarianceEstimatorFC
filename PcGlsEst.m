classdef PcGlsEst < handle
    properties
        N
        Nrf
        Q
        M
        Psi
        Sigma
        Pf
        Omega
        Ps
        Bm_dic
    end
    
    methods
        function obj = PcGlsEst(N, Nrf)
            % Constructor to initialize parameters and operators
            obj.N = N;
            obj.Nrf = Nrf;
            obj.Q = N/Nrf;
            obj.M = 2*obj.Q - 1;
            
            [obj.Psi, obj.Sigma, obj.Pf, obj.Omega, obj.Ps] = obj.create_linear_operators();
            [obj.Bm_dic] = obj.create_codebook();
        end
        
        function [Bm_dic] = create_codebook(obj)
            uu = 0 : obj.N-1;
            vv = 0 : obj.Nrf-1;
            Bm_dic = zeros(obj.N,obj.Nrf,obj.M);
            wm = 2*pi/(2*obj.Q-1)*(0:2*obj.Q-1);
            for m = 1 : obj.M
                for u = 1 : length(uu)
                    for v = 1 : length(vv)
                        if (vv(v)*obj.Q <= uu(u)) && ( uu(u) <= (vv(v)+1)*obj.Q-1) 
                            Bm_dic(u,v,m) = 1/sqrt(obj.Q)*exp(1j*wm(m)*(uu(u)-vv(v)*obj.Q));
                        end
                    end
                end
            end
        end
        
        function [Psi, Sigma, Pf, Omega, Ps] = create_linear_operators(obj)
            N_val = obj.N;
            Nrf_val = obj.Nrf;
            Q_val = obj.Q;
            M_val = obj.M;
            
            % --- Selection Matrix Ps ---
            Ps = zeros((2*Nrf_val-1)*(2*Q_val-1), 2*N_val-1);
            for q = 0 : 2*Nrf_val-2
                Ps(q*(2*Q_val-1)+1 : (q+1)*(2*Q_val-1), :) = [zeros(2*Q_val-1, q*Q_val), eye(2*Q_val-1), zeros(2*Q_val-1, 2*N_val-(q+2)*Q_val)];
            end
            
            % --- Spectral operator Omega ---
            Omega = zeros(2*Q_val-1, 2*Q_val-1);
            qq = -Q_val+1 : Q_val-1;
            uu = 0 : 2*Q_val-2;
            for u = 1 : length(uu)
                for q = 1 : length(qq)
                    Omega(u,q) = (1 - abs(qq(q))/Q_val) * exp(-1j*2*pi/(2*Q_val-1) * uu(u) * qq(q));
                end
            end
            
            % --- Spectral frequency selector Pf ---
            Pf = zeros(M_val*(2*Nrf_val-1), (2*Nrf_val-1)*(2*Q_val-1));
            for m = 1 : M_val
                em = zeros(2*Q_val-1, 1);
                em(m) = 1;
                Pf((m-1)*(2*Nrf_val-1)+1 : m*(2*Nrf_val-1), :) = kron(eye(2*Nrf_val-1), em.');
            end
            
            % --- Toeplitz structure operator Sigma ---
            Sigma = zeros(Nrf_val*Nrf_val, 2*Nrf_val-1);
            Sigma(:, Nrf_val) = reshape(eye(Nrf_val), [Nrf_val*Nrf_val, 1]);
            for i = 1 : Nrf_val-1
                Qi = [zeros(Nrf_val-i, i), eye(Nrf_val-i); zeros(i, Nrf_val-i), zeros(i, i)];
                Sigma(:, Nrf_val-i) = reshape(Qi, [Nrf_val*Nrf_val, 1]);
                Sigma(:, Nrf_val+i) = reshape(Qi.', [Nrf_val*Nrf_val, 1]); 
            end
            
            % --- Final Linear Operator Psi ---
            Psi = kron(eye(M_val), Sigma) * Pf * kron(eye(2*Nrf_val-1), Omega) * Ps;
        end
        
        function [hat_theta,hat_R] = estimate_doa(obj, x, L)
            % Y: Nrf x Km x M matrix
            % L: Number of sources to estimate
            
            [N_dim, K] = size(x);
            Km = K/obj.M;
            % Ensure the input Y matches the expected operator dimensions
            if  N_dim ~= obj.N
                error('Input matrix Y dimensions do not match initialized GlsEst parameters.');
            end

            x_reshaped = reshape(x,[obj.N,Km,obj.M]);
            inv_hat_Rep = zeros(obj.M*obj.Nrf*obj.Nrf,obj.M*obj.Nrf*obj.Nrf);
            p = zeros(obj.Nrf*obj.Nrf,obj.M);
            for m = 1 : obj.M
                ym = obj.Bm_dic(:,:,m)'*x_reshaped(:,:,m);
                hat_Sm = 1/Km*(ym*ym');
    
                
                inv_hat_Sm = eye(obj.Nrf)/hat_Sm;

                p(:,m) = hat_Sm(:);

                inv_hat_Repm = 1/Km*kron(inv_hat_Sm.',inv_hat_Sm);
                
                idx = (m-1)*obj.Nrf*obj.Nrf+1:m*obj.Nrf*obj.Nrf; 
                inv_hat_Rep(idx,idx) = inv_hat_Repm;
            end
            p = p(:);
            hat_r = (obj.Psi'*inv_hat_Rep*obj.Psi)\(obj.Psi'*inv_hat_Rep*p);
            hat_R = toeplitz(hat_r(obj.N:end),hat_r(obj.N:-1:1));
    
            hat_theta = rootmusic(hat_R,L);
 
        end
    end
end