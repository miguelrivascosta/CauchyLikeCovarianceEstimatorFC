

N = 8;
Nrf = [2,4];
% K = ceil(N/(Nrf-1))*100;
K = 240;
theta = deg2rad([-2.56,2.56]);
Rs = eye(length(theta));
snrdb = -20:2.5:10;
MC = 250;

run_rmse_vs_snr(N,Nrf(1),K,theta,Rs,snrdb,MC,method="clgls",rcrb="true")
run_rmse_vs_snr(N,Nrf(1),K,theta,Rs,snrdb,MC,method="bsa")
run_rmse_vs_snr(N,Nrf(1),K,theta,Rs,snrdb,MC,method="eaml")
run_rmse_vs_snr(N,Nrf(1),K,theta,Rs,snrdb,MC,method="ls2dft")
run_rmse_vs_snr(N,Nrf(1),K,theta,Rs,snrdb,MC,method="pi")
% run_rmse_vs_snr(N,Nrf,K,theta,Rs,snrdb,MC,method="dcomp")

run_rmse_vs_snr(N,Nrf(2),K,theta,Rs,snrdb,MC,method="clgls",rcrb="true")
run_rmse_vs_snr(N,Nrf(2),K,theta,Rs,snrdb,MC,method="bsa")
run_rmse_vs_snr(N,Nrf(2),K,theta,Rs,snrdb,MC,method="eaml")
run_rmse_vs_snr(N,Nrf(2),K,theta,Rs,snrdb,MC,method="ls2dft")
run_rmse_vs_snr(N,Nrf(2),K,theta,Rs,snrdb,MC,method="pi")


