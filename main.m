clear, close all

N = 256;
Nrf = 64;

clgls = CLGLS(N,Nrf);

Km = 100;
K = clgls.M*Km;

theta_deg = [-30,45];
theta = deg2rad(theta_deg);
L = length(theta);

A = exp(1j*pi*(0:N-1).'*sin(theta));
Rs = eye(L);
snrdb = -20:1:10;

MC = 500;

rmse   = zeros(1, length(snrdb));

for i = 1 : length(snrdb)
    sigman2 = 10.^(-snrdb(i)/10);
    R = A*Rs*A'+sigman2*eye(N);
    sum_se = 0;
    for mc = 1 : MC
        x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,sigman2,K);
        hat_theta = clgls.estimate_doas(x,L,mode="fast");
        err   = theta - hat_theta;
        sum_se = sum_se + mean(err.^2);
    end
    rmse(i)   = rad2deg(sqrt(sum_se / MC));
    fprintf("Iteration: %d/%d; RMSE=%f\n",i,length(snrdb),rmse(i))

end

snrdb_rcrb = snrdb(1):0.1:snrdb(end);
obj = CLGLS(N,Nrf);
crb = ToolsObj.eval_crb_theta_per_snr(1/2,1,theta,power(10,snrdb_rcrb/10),obj.Bm_dic,Km);
rcrb_deg = rad2deg(sqrt(mean(crb,1)));

figure
hold on
semilogy(snrdb,rmse)
semilogy(snrdb_rcrb,rcrb_deg)
hold off
set(gca,yscale='log')

