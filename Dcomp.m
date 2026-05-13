classdef Dcomp < handle
    properties
        N           % Number of BS antennas
        Nrf         % Number of RF chains
        Q
        T           % Number of snapshots / beamforming matrices
        D           % Dictionary size (grid size for AoA)
        theta_grid_deg
        theta_grid  % Spatial grid
        A_dict      % Dictionary matrix
        Phi_cell    % Effective sensing matrices (Psi_t)
        Bm_dic    % Beamforming matrices (W_t)
        Phi_dic
        M
    end
    
    methods
        function obj = Dcomp(N, Nrf, T, D)
            % CONSTRUCTOR: Initializes parameters and builds the codebook/dictionary
            obj.N = N;
            obj.Nrf = Nrf;
            obj.Q = obj.N/obj.Nrf;
            obj.T = T;   % Note: Represents 'M' (beamforming matrices) from your prompt
            obj.M = T;
            obj.D = D;
            
            % 1. Create Spatial Grid and Dictionary Matrix
            obj.theta_grid_deg = linspace(-90, 90, D);
            obj.theta_grid = deg2rad(obj.theta_grid_deg);

            obj.A_dict = zeros(N, D);
            for d = 1:D
                obj.A_dict(:, d) = exp(1j * pi * (0:N-1)' .* sin(obj.theta_grid(d)));
            end
            
            % 2. Create Codebook and Effective Sensing Matrices (Phi_cell)
            Q = N / Nrf;
            obj.Phi_cell = cell(1, T);
            obj.Bm_dic = zeros(obj.N, obj.Nrf, T);
            
            for t = 1:T
                % A. Time-Varying Analog Combiner (Random Phase Shifters)
                W_t = exp(1j*unifrnd(0,2*pi, [obj.N, obj.Nrf])) / sqrt(obj.N);
                
                % B. Store Codebook and Calculate Effective Sensing Matrix
                obj.Bm_dic(:,:,t) = W_t;
                obj.Phi_cell{t} = W_t' * obj.A_dict;

                obj.Phi_dic(:,:,t) = W_t' * obj.A_dict;
            end
        end
        
        function update_dictionaries(obj)
            obj.Phi_cell = cell(1, obj.T);
            obj.Bm_dic = zeros(obj.N, obj.Nrf, obj.T);
            
            for t = 1:obj.T
                % A. Time-Varying Analog Combiner (Random Phase Shifters)
                W_t = exp(1j*unifrnd(0,2*pi, [obj.N, obj.Nrf])) / sqrt(obj.N);

                % B. Store Codebook and Calculate Effective Sensing Matrix
                obj.Bm_dic(:,:,t) = W_t;
                obj.Phi_cell{t} = W_t' * obj.A_dict;

                obj.Phi_dic(:,:,t) = W_t' * obj.A_dict;
            end
        end
        function [hat_theta,Rg_S] = estimate_doas(obj, x, L)
            % Pre-allocate 3D arrays for vectorized access
            % Ry = zeros(obj.Nrf, obj.Nrf, obj.T);
            % for t = 1:obj.T
            %     Y_t = obj.Bm_dic(:,:,t)' * x(:,t);
            %     Ry(:,:,t) = Y_t * Y_t';
            % end
            Y = pagemtimes(pagectranspose(obj.Bm_dic),reshape(x,[obj.N,1,size(x,2)]));
            Ry = pagemtimes(Y,pagectranspose(Y));
            V = Ry; % Initialize Residuals
            S = []; 
            
            % OMP Loop
            for l = 1:L
                % VECTORIZED METRIC CALCULATION
                % Instead of looping over D, we compute the diagonal of (Phi' * V * Phi)
                % across all snapshots t simultaneously.
                % metric = zeros(1, obj.D);

                VPhi = pagemtimes(V,obj.Phi_dic);
                metric = real(sum(conj(obj.Phi_dic) .* VPhi, [1, 3]));
                % for t = 1:obj.T
                %     % VPhi results in (Nrf x D)
                %     % VPhi = V(:,:,t) * obj.Phi_cell{t};
                %     % The quadratic form phi' * V * phi is the column-wise dot product
                %     metric = metric + real(sum(conj(obj.Phi_dic(:,:,t)) .* VPhi(:,:,t), 1));
                % end
                
                if ~isempty(S), metric(S) = -inf; end
                
                [~, j] = max(metric);
                S = [S, j];
                
                % FASTER RESIDUAL UPDATE
                % % Using QR decomposition to find the orthogonal complement projection
                % for t = 1:obj.T
                %     Phi_t_S = obj.Phi_cell{t}(:, S);
                %     [Q_orth, ~] = qr(Phi_t_S, 0); 
                %     % Orthogonal projection: P_perp = I - Q*Q'
                %     % V = P_perp * Ry * P_perp
                %     % Optimized as: Ry_proj = Ry - Q*(Q'*Ry) -> then same for the right side
                %     temp = Ry(:,:,t) - Q_orth * (Q_orth' * Ry(:,:,t));
                %     V(:,:,t) = temp - (temp * Q_orth) * Q_orth';
                % end
                if ~(size(S,2) == L)
                    for t = 1:obj.T
                        Phi_t_S = obj.Phi_cell{t}(:, S);
                        [Q_orth, ~] = qr(Phi_t_S, 0);  % Phi_t_S = Q·R, so P = Q·Q^H
                    
                        % CORRECT: Ry - P·Ry·P  (Algorithm 5, COMP/DCOMP formula)
                        QHRyQ = Q_orth' * Ry(:,:,t) * Q_orth;      % Q^H · Ry · Q  (small: |S|×|S|)
                        V(:,:,t) = Ry(:,:,t) - Q_orth * QHRyQ * Q_orth';  % Ry - Q·(Q^H·Ry·Q)·Q^H
                    end
                end
            end
            
            % FINAL EXTRACTION
            Rg_S = zeros(obj.D, obj.D); 
            for t = 1:obj.T
                Phi_t_S = obj.Phi_cell{t}(:, S);
                % Use backslash (\) instead of pinv for speed and stability
                % pinv(A)*B is roughly A\B
                sol = Phi_t_S \ Ry(:,:,t);
                Rg_S(S,S) = Rg_S(S,S) + (sol / (Phi_t_S')); 
            end
            
            [~, idx] = maxk(real(diag(Rg_S)), L);
            hat_theta = obj.theta_grid(sort(idx, 'ascend'));

            % hat_Rx = obj.A_dict*Rg_S*obj.A_dict';
            % hat_theta = sort(rootmusic(hat_Rx,L),'ascend');
        end
    end
end