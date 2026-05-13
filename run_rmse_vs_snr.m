function run_rmse_vs_snr(N,Nrf,K,theta,Rs,snrdb,MC,options)
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
        obj = SheinvaldMethod(N,Nrf,-90,90,.1);
    elseif options.method == "pi"
        fprintf("\n Simulating RMSE vs SNR with Phase-Independent\n")
        obj = PhaseIndependent(N,Nrf);
    elseif options.method == "bsa"
        fprintf("\n Simulating RMSE vs SNR with BSA\n")
        obj = BSA(N,Nrf);
    else
        fprintf("\n Simulating RMSE vs SNR with DCOMP\n")
        obj = Dcomp(N,Nrf,K,1000);
    end
    
    filename = sprintf('rmse_vs_snr_%s_n%d_Nrf%d.mat', options.method, N, Nrf);

    if options.rcrb == "true"
        snrdb_rcrb = snrdb(1):0.1:snrdb(end);
        crb = ToolsObj.eval_crb_theta_per_snr(1/2,1,theta,power(10,snrdb_rcrb/10),obj.Bm_dic,K/obj.M);
        rcrb_deg = rad2deg(sqrt(mean(crb,1)));
    else
        snrdb_rcrb = NaN;
        rcrb_deg = NaN;
    end



    rmse_deg = zeros(size(snrdb));

    A = exp(1j*pi*(0:N-1).'*sin(theta));
    L = length(theta);
    for i = 1 : length(snrdb)
        sigman2 = 10.^(-snrdb(i)/10);
        sum_se = 0;

        for mc = 1 : MC
            if options.method == "dcomp"
                obj.update_dictionaries();
            end
            x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,sigman2,K);
            hat_theta = obj.estimate_doas(x,L);
            err   = theta - hat_theta;
            sum_se = sum_se + mean(err.^2);
        end

        rmse_deg(i)   = rad2deg(sqrt(sum_se / MC));
        fprintf("Iteration: %d/%d; RMSE=%f\n",i,length(snrdb),rmse_deg(i))
    end

    data = struct(...
        'N', N, ...
        'Nrf', Nrf, ...
        'Km', K/obj.M, ...
        'M', obj.M, ...
        'K', K, ...
        'theta', theta, ...
        'rmse_deg', rmse_deg, ...
        'snrdb', snrdb, ...
        'rcrb_deg', rcrb_deg, ...
        'snrdb_rcrb', snrdb_rcrb, ...
        'MC', MC);

    if isfile(filename)
        filename = "new_" + filename;
    end
    
    %semilogy(data.snrdb,data.rmse)
    save(filename, '-struct','data');

end