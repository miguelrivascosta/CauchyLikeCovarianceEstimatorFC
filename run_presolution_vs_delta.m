function run_presolution_vs_delta(N,Nrf,K,delta_th,Rs,snrdb,MC,options)
    arguments
        N
        Nrf
        K
        delta_th
        Rs
        snrdb
        MC
        options.method (1,1) string {mustBeMember(options.method, ["clgls", "ls2dft","eaml", "pi","bsa","dcomp"])}
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
    
    fprintf("Simulating Presolution vs delta. N=%d, Nrf=%d, Method->%s\n",N,Nrf,options.method)

    filename = sprintf('presolution_vs_delta_%s_n%d_Nrf%d.mat', options.method, N, Nrf);
    

    presolution = zeros(size(delta_th));
   
    for i = 1 : length(delta_th)
        sum_presolution = 0;
        theta = [0,delta_th(i)];
        A = exp(1j*pi*(0:N-1).'*sin(theta));

        for mc = 1 : MC
            if options.method == "dcomp"
                obj.create_dictionaries(K);
            end
            x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,10.^(-snrdb/10),K);
            hat_theta = obj.estimate_doas(x,2);
            
            mid_point = (abs(theta(1)-theta(2))/2);
            test_theta1 = (abs(hat_theta(1) - theta(1))) < mid_point;
            test_theta2 = (abs(hat_theta(2) - theta(2))) < mid_point;
            if test_theta1 && test_theta2
                sum_presolution = sum_presolution+1;
            end
            
        end

        presolution(i)   = sum_presolution / MC;

        if options.show_iter == "true"
            fprintf("Iteration: %d/%d; RMSE=%f\n",i,length(delta_th),presolution(i))
        end
    end

    data = struct(...
        'N', N, ...
        'Nrf', Nrf, ...
        'Km', K/obj.M, ...
        'M', obj.M, ...
        'K', K, ...
        'presolution', presolution, ...
        'delta_deg', rad2deg(delta_th), ...
        'MC', MC);

    if isfile(filename)
        filename = "new_" + filename;
    end
    
    %semilogy(data.snrdb,data.rmse)
    save(filename, '-struct','data');

end