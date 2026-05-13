clear, close all


color_clgls="blue";
line_clgls = "-";
name_clgls="ClGLS";

color_bsa="red";
line_bsa = "-";
name_bsa="BSA";

color_pi="green";
line_pi = "-";
marker_pi = "o";
name_pi="PI";

color_eaml="yellow";
line_eaml = "-";
name_eaml="EAML";

color_ls2dft="magenta";
line_ls2dft = "-";
name_ls2dft="LS-2DFT";

color_dcomp="cyan";
line_dcomp = "-";
name_dcomp="DCOMP";

color_rcrb="black";
line_rcrb = "-";
marker_rcrb = "none";
name_rcrb="RCRB";




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% RMSE vs SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dcomp_Nrf2 = load('rmse_vs_snr_dcomp_n8_Nrf2.mat');
ls2dft_Nrf2 = load('rmse_vs_snr_ls2dft_n8_Nrf2.mat');
eaml_Nrf2 = load('rmse_vs_snr_eaml_n8_Nrf2.mat');
pi_Nrf2 = load('rmse_vs_snr_pi_n8_Nrf2.mat');
bsa_Nrf2 = load('rmse_vs_snr_bsa_n8_Nrf2.mat');
clgls_Nrf2 = load('rmse_vs_snr_clgls_n8_Nrf2.mat');

% dcomp_Nrf4 = load('rmse_vs_snr_dcomp_n8_Nrf4.mat');
ls2dft_Nrf4 = load('rmse_vs_snr_ls2dft_n8_Nrf4.mat');
eaml_Nrf4 = load('rmse_vs_snr_eaml_n8_Nrf4.mat');
pi_Nrf4 = load('rmse_vs_snr_pi_n8_Nrf4.mat');
bsa_Nrf4 = load('rmse_vs_snr_bsa_n8_Nrf4.mat');
clgls_Nrf4 = load('rmse_vs_snr_clgls_n8_Nrf4.mat');

marker_Nrf2 = "o";
marker_Nrf4 = ">";

figure
hold on
% semilogy(dcomp_Nrf2.snrdb,dcomp_Nrf2.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp)
semilogy(ls2dft_Nrf2.snrdb,ls2dft_Nrf2.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft)
semilogy(eaml_Nrf2.snrdb,eaml_Nrf2.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml)
semilogy(pi_Nrf2.snrdb,pi_Nrf2.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi)
semilogy(bsa_Nrf2.snrdb,bsa_Nrf2.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa)
semilogy(clgls_Nrf2.snrdb,clgls_Nrf2.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls)
semilogy(clgls_Nrf2.snrdb_rcrb,clgls_Nrf2.rcrb_deg,LineStyle=line_rcrb,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb)

% semilogy(dcomp_Nrf4.snrdb,dcomp_Nrf4.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp)
semilogy(ls2dft_Nrf4.snrdb,ls2dft_Nrf4.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft)
semilogy(eaml_Nrf4.snrdb,eaml_Nrf4.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml)
semilogy(pi_Nrf4.snrdb,pi_Nrf4.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi)
semilogy(bsa_Nrf4.snrdb,bsa_Nrf4.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa)
semilogy(clgls_Nrf4.snrdb,clgls_Nrf4.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls)
semilogy(clgls_Nrf4.snrdb_rcrb,clgls_Nrf4.rcrb_deg,LineStyle=line_rcrb,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb)

hold off
set(gca,yscale='log')
xlabel('$\mathrm{SNR [dB]}$',Interpreter='latex')
ylabel('$\mathrm{RMSE [deg]}$',Interpreter='latex')
legend(Interpreter='latex')


