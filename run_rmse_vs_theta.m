function run_rmse_vs_theta(N,Nrf,K,theta,Rs,snrdb,MC,options)
    arguments
        N
        Nrf
        K
        theta
        Rs
        snrdb
        MC
        options.method (1,1) string {mustBeMember(options.method, ["clgls", "ls2dft","eaml", "pi","bsa","dcomp"])}
        options.rcrb (1,1) string {mustBeMember(options.rcrb, ["true", "false"])} = "false"
    end
    
    if options.method == "clgls"
        fprintf("\n Simulating RMSE vs SNR with ClGLS\n")
        obj = CLGLS(N,Nrf);
    elseif options.method == "ls2dft"
        fprintf("\n Simulating RMSE vs SNR with LS-2DFT\n")
        obj = LS_2DFT(N,Nrf);
    elseif options.method == "eaml"
        fprintf("\n Simulating RMSE vs SNR with E-AML\n")
        obj = EAML(N,Nrf);
    elseif options.method == "pi"
        fprintf("\n Simulating RMSE vs SNR with Phase-Independent\n")
        obj = PhaseIndependent(N,Nrf);
    elseif options.method == "bsa"
        fprintf("\n Simulating RMSE vs SNR with BSA\n")
        obj = BSA(N,Nrf);
    else
        fprintf("\n Simulating RMSE vs SNR with DCOMP\n")
        obj = Dcomp(N,Nrf);
    end
    
    filename = sprintf('rmse_vs_theta_%s_n%d_Nrf%d.mat', options.method, N, Nrf);
    
    sigman2 = 10.^(-snrdb/10);

    if options.rcrb == "true"
        theta_rcrb_deg = theta(1):deg2rad(0.1):theta(end);
        rcrb_deg = zeros(size(delta_rcrb));
        for i = 1 : length(delta_rcrb)
            aux = ToolsObj.crb_unc_uncorrelated_theta(1/2,1,theta(i),Rs,sigman2,obj.Bm_dic,K/obj.M);
            crb = diag(aux(1:2,1:2));
            rcrb_deg(i) = rad2deg(sqrt(mean(crb,1)));
        end
    else
        theta_rcrb_deg = NaN;
        rcrb_deg = NaN;
    end

    rmse_deg = zeros(size(theta));
    L = 1;

   

    for i = 1 : length(theta)
        sum_se = 0;
        A = exp(1j*pi*(0:N-1).'*sin(theta(i)));
        for mc = 1 : MC
            if options.method == "dcomp"
                obj.create_dictionaries(K);
            end
            x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,sigman2,K);
            hat_theta = obj.estimate_doas(x,L);
            err   = theta(i) - hat_theta;
            sum_se = sum_se + mean(err.^2);
        end

        rmse_deg(i)   = rad2deg(sqrt(sum_se / MC));
        fprintf("Iteration: %d/%d; RMSE=%f\n",i,length(theta),rmse_deg(i))
    end

    data = struct(...
        'N', N, ...
        'Nrf', Nrf, ...
        'Km', K/obj.M, ...
        'M', obj.M, ...
        'K', K, ...
        'rmse_deg', rmse_deg, ...
        'theta_deg', rad2deg(theta), ...
        'rcrb_deg', rcrb_deg, ...
        'theta_rcrb_deg', rad2deg(theta_rcrb_deg), ...
        'MC', MC);

    if isfile(filename)
        filename = "new_" + filename;
    end
    
    %semilogy(data.theta_deg,data.rmse_deg)
    save(filename, '-struct','data');

end