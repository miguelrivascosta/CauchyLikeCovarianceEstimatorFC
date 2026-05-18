function run_rmse_vs_K(N,Nrf,K,theta,Rs,snrdb,MC,options)
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
        fprintf("\n Simulating RMSE vs K with ClGLS\n")
        obj = CLGLS(N,Nrf);
    elseif options.method == "ls2dft"
        fprintf("\n Simulating RMSE vs K with LS-2DFT\n")
        obj = LS_2DFT(N,Nrf);
    elseif options.method == "eaml"
        fprintf("\n Simulating RMSE vs K with E-AML\n")
        obj = EAML(N,Nrf);
    elseif options.method == "pi"
        fprintf("\n Simulating RMSE vs K with Phase-Independent\n")
        obj = PhaseIndependent(N,Nrf);
    elseif options.method == "bsa"
        fprintf("\n Simulating RMSE vs K with BSA\n")
        obj = BSA(N,Nrf);
    else
        fprintf("\n Simulating RMSE vs K with DCOMP\n")
        obj = Dcomp(N,Nrf);
    end
    
    filename = sprintf('rmse_vs_K_%s_n%d_Nrf%d.mat', options.method, N, Nrf);
    
    sigman2 = 10.^(-snrdb/10);

    if options.rcrb == "true"
        K_rcrb = K;
        rcrb_deg = zeros(size(K_rcrb));
        for i = 1 : length(K_rcrb)
            aux = ToolsObj.crb_unc_uncorrelated_theta(1/2,1,theta,Rs,sigman2,obj.Bm_dic,K_rcrb(i)/obj.M);
            crb = diag(aux(1:2,1:2));
            rcrb_deg(i) = rad2deg(sqrt(mean(crb,1)));
        end
    else
        K_rcrb = NaN;
        rcrb_deg = NaN;
    end

    rmse_deg = zeros(size(K));
    L = 2;

   

    for i = 1 : length(K)
        sum_se = 0;
        A = exp(1j*pi*(0:N-1).'*sin(theta));
        for mc = 1 : MC
            if options.method == "dcomp"
                obj.create_dictionaries(K(i));
            end
            x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,sigman2,K(i));
            hat_theta = obj.estimate_doas(x,L);
            err   = theta - hat_theta;
            sum_se = sum_se + mean(err.^2);
        end

        rmse_deg(i)   = rad2deg(sqrt(sum_se / MC));
        fprintf("Iteration: %d/%d; RMSE=%f\n",i,length(K),rmse_deg(i))
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
        'K_rcrb', K_rcrb, ...
        'MC', MC);

    if isfile(filename)
        filename = "new_" + filename;
    end
    
    %semilogy(data.snrdb,data.rmse)
    save(filename, '-struct','data');

end