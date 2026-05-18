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
        options.show_iter (1,1) string {mustBeMember(options.show_iter, ["true", "false"])} = "false"
    end
    
    if options.method == "clgls"
        obj = CLGLS(N,Nrf);
    elseif options.method == "ls2dft"
        obj = LS_2DFT(N,Nrf);
    elseif options.method == "eaml"
        obj = EAML(N,Nrf);
    elseif options.method == "pi"
        obj = PhaseIndependent(N,Nrf);
    elseif options.method == "bsa"
        obj = BSA(N,Nrf);
    else
        obj = Dcomp(N,Nrf);
    end
    
    fprintf("Simulating RMSE vs SNR. N=%d, Nrf=%d, Method->%s\n",N,Nrf,options.method)

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
                obj.create_dictionaries(K);
            end
            x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,sigman2,K);
            hat_theta = obj.estimate_doas(x,L);
            err   = theta - hat_theta;
            sum_se = sum_se + mean(err.^2);
        end

        rmse_deg(i)   = rad2deg(sqrt(sum_se / MC));
        if options.show_iter == "true"
            fprintf("Iteration: %d/%d; RMSE=%f\n",i,length(snrdb),rmse_deg(i))
        end
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
    
    %semilogy(data.snrdb,data.rmse_deg)
    save(filename, '-struct','data');

end