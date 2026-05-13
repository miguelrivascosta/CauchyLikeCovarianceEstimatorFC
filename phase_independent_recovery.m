function [Ry_hat] = phase_independent_recovery(Y_dig,Nrf,K)
    %%% Y: N x T matrix where each column is the received signal y = As+n
    %%% K: 2K is the min number of beams to recover the covariance
    %%% Q: Number of snapshots to average each beam
    %%% Reference: 
    % O. Shmonin et al., "Phase-Independent Beamspace MUSIC Algorithm for a Single Port Phased Antenna Array," 
    % 2023 31st Telecommunications Forum (TELFOR), Belgrade, Serbia, 2023, pp. 1-4, doi: 10.1109/TELFOR59449.2023.10372794.
    
    [N,size_data] = size(Y_dig);
    M = 2*K;

    num_rep = size_data/(M/Nrf);
    k = -K:K-1;
    psi_k = k*pi/K;
    B = 1/sqrt(N)*ToolsObj.a_ula(N,psi_k);
    
    codebook_size = M/Nrf;
    codebook = zeros(N,Nrf,codebook_size);
    idx_mat = zeros(Nrf,codebook_size);
    % idx = 1:Nrf;
    for i = 1 : codebook_size
        idx = i : 2*N/Nrf : 2*N;
        idx_mat(:,i) = idx.';
        codebook(:,:,i) = B(:,idx);
        % idx = idx_mat(:,i) + Nrf;
    end
    % aux = reshape(Y_dig,[N,Q,M/Nrf]);


    avg_abs_y = zeros(2*N,1);
    t = 1;
    for i = 1 : codebook_size
        Bm = codebook(:,:,i);
        for j = 1 : num_rep
            idx = idx_mat(:,i);
            avg_abs_y(idx) = avg_abs_y(idx) + power(abs(Bm'*Y_dig(:,t)),2);
            t = t + 1;
        end
    end
    avg_abs_y = avg_abs_y / num_rep;
    % Y_hyb = pagemtimes(pagectranspose(B),aux);
    % avg_abs_y = N*mean(power(abs(Y_hyb),2),2);
    % abs_y = reshape(mean(power(abs(Y_hyb),2),2),[M,1]);
    % stem(abs_y)
    f_hat = reshape(avg_abs_y,[M,1]);

    q = 1-N:N-1;
    exp_fourier_series = exp(-1j*pi*k.'*q/K);
    Dq = zeros(2*N-1,1);    
    rq = zeros(length(q),1);
    for ii = 1:length(q)
        Dq(ii) = 1/(2*K)*sum(f_hat.*exp_fourier_series(:,ii),1);
        rq(ii) = Dq(ii)/(N-abs(q(ii)));
    end
    c = rq(N:-1:1);
    r = rq(N:1:end);
    Ry_hat = toeplitz(c,r);
end