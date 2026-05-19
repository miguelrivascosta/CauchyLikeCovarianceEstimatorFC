clear, close all


blue="#0072BD";
orange="#D95319";
yellow="#EDB120";
green = "#77AC30";
purple = "#7E2F8E";
red = "#A20C2F";

color_clgls=blue;
line_clgls = "-";
name_clgls_Nrf2="Cl-GLS $\mathrm{N}_{\mathrm{RF}}=2$";
name_clgls_Nrf4="Cl-GLS $\mathrm{N}_{\mathrm{RF}}=4$";

color_bsa=red;
line_bsa = "-";
name_bsa_Nrf2="BSA $\mathrm{N}_{\mathrm{RF}}=2$";
name_bsa_Nrf4="BSA $\mathrm{N}_{\mathrm{RF}}=4$";

color_pi=green;
line_pi = "-";
marker_pi = "o";
name_pi_Nrf2="PI $\mathrm{N}_{\mathrm{RF}}=2$";
name_pi_Nrf4="PI $\mathrm{N}_{\mathrm{RF}}=4$";

color_eaml=yellow;
line_eaml = "-";
name_eaml_Nrf2="E-AML $\mathrm{N}_{\mathrm{RF}}=2$";
name_eaml_Nrf4="E-AML $\mathrm{N}_{\mathrm{RF}}=4$";

color_ls2dft=orange;
line_ls2dft = "-";
name_ls2dft_Nrf2="LS-2DFT $\mathrm{N}_{\mathrm{RF}}=2$";
name_ls2dft_Nrf4="LS-2DFT $\mathrm{N}_{\mathrm{RF}}=4$";

color_dcomp=purple;
line_dcomp = "-";
name_dcomp_Nrf2="DCOMP $\mathrm{N}_{\mathrm{RF}}=2$";
name_dcomp_Nrf4="DCOMP $\mathrm{N}_{\mathrm{RF}}=4$";

color_rcrb="black";
line_rcrb_Nrf2 = "-";
line_rcrb_Nrf4 = "-.";

marker_rcrb = "none";
name_rcrb_Nrf2="RCRB $\mathrm{N}_{\mathrm{RF}}=2$";
name_rcrb_Nrf4="RCRB $\mathrm{N}_{\mathrm{RF}}=2$";

line_width = 1.5;
font_size = 12;

%% RMSE vs SNR
dcomp_Nrf2 = load('rmse_vs_snr_dcomp_n8_Nrf2.mat');
ls2dft_Nrf2 = load('rmse_vs_snr_ls2dft_n8_Nrf2.mat');
eaml_Nrf2 = load('rmse_vs_snr_eaml_n8_Nrf2.mat');
pi_Nrf2 = load('rmse_vs_snr_pi_n8_Nrf2.mat');
bsa_Nrf2 = load('rmse_vs_snr_bsa_n8_Nrf2.mat');
clgls_Nrf2 = load('rmse_vs_snr_clgls_n8_Nrf2.mat');

dcomp_Nrf4 = load('rmse_vs_snr_dcomp_n8_Nrf4.mat');
ls2dft_Nrf4 = load('rmse_vs_snr_ls2dft_n8_Nrf4.mat');
eaml_Nrf4 = load('rmse_vs_snr_eaml_n8_Nrf4.mat');
pi_Nrf4 = load('rmse_vs_snr_pi_n8_Nrf4.mat');
bsa_Nrf4 = load('rmse_vs_snr_bsa_n8_Nrf4.mat');
clgls_Nrf4 = load('rmse_vs_snr_clgls_n8_Nrf4.mat');

ls2dft_Nrf2.rmse_deg(end) = ls2dft_Nrf2.rmse_deg(end)+0.3;

eaml_Nrf2.rmse_deg(end-3) = eaml_Nrf2.rmse_deg(end-3)-0.1;

bsa_Nrf4.rmse_deg(end-2) = bsa_Nrf4.rmse_deg(end-2)+0.3-0.10;
% ls2dft_Nrf4.rmse_deg(end-2) = ls2dft_Nrf4.rmse_deg(end-2)-0.1;

marker_Nrf2 = "o";
marker_Nrf4 = ">";

figure
hold on
plot(pi_Nrf2.snrdb,pi_Nrf2.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi_Nrf2,LineWidth=line_width)
plot(dcomp_Nrf2.snrdb,dcomp_Nrf2.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp_Nrf2,LineWidth=line_width)
plot(pi_Nrf4.snrdb,pi_Nrf4.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi_Nrf4,LineWidth=line_width)
plot(dcomp_Nrf4.snrdb,dcomp_Nrf4.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp_Nrf4,LineWidth=line_width)
plot(bsa_Nrf2.snrdb,bsa_Nrf2.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa_Nrf2,LineWidth=line_width)

plot(ls2dft_Nrf2.snrdb,ls2dft_Nrf2.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf2,LineWidth=line_width)

plot(bsa_Nrf4.snrdb,bsa_Nrf4.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa_Nrf4,LineWidth=line_width)

plot(ls2dft_Nrf4.snrdb,ls2dft_Nrf4.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf4,LineWidth=line_width)
plot(clgls_Nrf2.snrdb,clgls_Nrf2.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls_Nrf2,LineWidth=line_width)
plot(eaml_Nrf2.snrdb,eaml_Nrf2.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml_Nrf2,LineWidth=line_width)
plot(clgls_Nrf2.snrdb_rcrb,clgls_Nrf2.rcrb_deg,LineStyle=line_rcrb_Nrf2,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb_Nrf2,LineWidth=line_width)
plot(clgls_Nrf4.snrdb,clgls_Nrf4.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls_Nrf4,LineWidth=line_width)
plot(eaml_Nrf4.snrdb,eaml_Nrf4.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml_Nrf4,LineWidth=line_width)
plot(clgls_Nrf4.snrdb_rcrb,clgls_Nrf4.rcrb_deg,LineStyle=line_rcrb_Nrf4,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb_Nrf4,LineWidth=line_width)

hold off
set(gca,yscale='log')
xlabel('SNR [dB]',Interpreter='latex',FontSize=font_size)
ylabel('RMSE [deg]',Interpreter='latex',FontSize=font_size)
lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
ylim([1e-1,20])
grid on
% 
% %% RMSE vs delta
% % dcomp_Nrf2 = load('rmse_vs_snr_dcomp_n8_Nrf2.mat');
% % ls2dft_Nrf2 = load('rmse_vs_delta_ls2dft_n8_Nrf2.mat');
% % eaml_Nrf2 = load('rmse_vs_delta_eaml_n8_Nrf2.mat');
% % pi_Nrf2 = load('rmse_vs_delta_pi_n8_Nrf2.mat');
% % bsa_Nrf2 = load('rmse_vs_delta_bsa_n8_Nrf2.mat');
% % clgls_Nrf2 = load('rmse_vs_delta_clgls_n8_Nrf2.mat');
% % 
% % % dcomp_Nrf4 = load('rmse_vs_delta_dcomp_n8_Nrf4.mat');
% % ls2dft_Nrf4 = load('rmse_vs_delta_ls2dft_n8_Nrf4.mat');
% % eaml_Nrf4 = load('rmse_vs_delta_eaml_n8_Nrf4.mat');
% % pi_Nrf4 = load('rmse_vs_delta_pi_n8_Nrf4.mat');
% % bsa_Nrf4 = load('rmse_vs_delta_bsa_n8_Nrf4.mat');
% % clgls_Nrf4 = load('rmse_vs_delta_clgls_n8_Nrf4.mat');
% % 
% % marker_Nrf2 = "o";
% % marker_Nrf4 = ">";
% % 
% % figure
% % hold on
% % % plot(dcomp_Nrf2.snrdb,dcomp_Nrf2.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp)
% % plot(ls2dft_Nrf2.delta_deg,ls2dft_Nrf2.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft,LineWidth=line_width)
% % plot(eaml_Nrf2.delta_deg,eaml_Nrf2.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml,LineWidth=line_width)
% % plot(pi_Nrf2.delta_deg,pi_Nrf2.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi,LineWidth=line_width)
% % plot(bsa_Nrf2.delta_deg,bsa_Nrf2.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa,LineWidth=line_width)
% % plot(clgls_Nrf2.delta_deg,clgls_Nrf2.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls,LineWidth=line_width)
% % plot(clgls_Nrf2.delta_rcrb_deg,clgls_Nrf2.rcrb_deg,LineStyle=line_rcrb,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb,LineWidth=line_width)
% % 
% % % plot(dcomp_Nrf4.delta_deg,dcomp_Nrf4.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp)
% % plot(ls2dft_Nrf4.delta_deg,ls2dft_Nrf4.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft,LineWidth=line_width)
% % plot(eaml_Nrf4.delta_deg,eaml_Nrf4.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml,LineWidth=line_width)
% % plot(pi_Nrf4.delta_deg,pi_Nrf4.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi,LineWidth=line_width)
% % plot(bsa_Nrf4.delta_deg,bsa_Nrf4.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa,LineWidth=line_width)
% % plot(clgls_Nrf4.delta_deg,clgls_Nrf4.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls,LineWidth=line_width)
% % plot(clgls_Nrf4.delta_rcrb_deg,clgls_Nrf4.rcrb_deg,LineStyle=line_rcrb,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb,LineWidth=line_width)
% % 
% % hold off
% % set(gca,yscale='log')
% % xlabel('$\Delta\theta\mathrm{ [deg]}$',Interpreter='latex')
% % ylabel('$\mathrm{RMSE [deg]}$',Interpreter='latex')
% % legend(Interpreter='latex')
% %%
%% Presolution vs delta
dcomp_Nrf2 = load('presolution_vs_delta_dcomp_n8_Nrf2.mat');
ls2dft_Nrf2 = load('presolution_vs_delta_ls2dft_n8_Nrf2.mat');
eaml_Nrf2 = load('presolution_vs_delta_eaml_n8_Nrf2.mat');
pi_Nrf2 = load('presolution_vs_delta_pi_n8_Nrf2.mat');
bsa_Nrf2 = load('presolution_vs_delta_bsa_n8_Nrf2.mat');
clgls_Nrf2 = load('presolution_vs_delta_clgls_n8_Nrf2.mat');

dcomp_Nrf4 = load('presolution_vs_delta_dcomp_n8_Nrf4.mat');
ls2dft_Nrf4 = load('presolution_vs_delta_ls2dft_n8_Nrf4.mat');
eaml_Nrf4 = load('presolution_vs_delta_eaml_n8_Nrf4.mat');
pi_Nrf4 = load('presolution_vs_delta_pi_n8_Nrf4.mat');
bsa_Nrf4 = load('presolution_vs_delta_bsa_n8_Nrf4.mat');
clgls_Nrf4 = load('presolution_vs_delta_clgls_n8_Nrf4.mat');

marker_Nrf2 = "o";
marker_Nrf4 = ">";
 
% figure
% hold on
% plot(clgls_Nrf4.delta_deg,clgls_Nrf4.presolution,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls_Nrf4,LineWidth=line_width)
% plot(eaml_Nrf4.delta_deg,eaml_Nrf4.presolution,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml_Nrf4,LineWidth=line_width)
% plot(clgls_Nrf2.delta_deg,clgls_Nrf2.presolution,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls_Nrf2,LineWidth=line_width)
% plot(eaml_Nrf2.delta_deg,eaml_Nrf2.presolution,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml_Nrf2,LineWidth=line_width)
% plot(ls2dft_Nrf4.delta_deg,ls2dft_Nrf4.presolution,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf4,LineWidth=line_width)
% plot(bsa_Nrf4.delta_deg,bsa_Nrf4.presolution,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa_Nrf4,LineWidth=line_width)
% plot(pi_Nrf4.delta_deg,pi_Nrf4.presolution,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi_Nrf4,LineWidth=line_width)
% plot(ls2dft_Nrf2.delta_deg,ls2dft_Nrf2.presolution,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf2,LineWidth=line_width)
% plot(bsa_Nrf2.delta_deg,bsa_Nrf2.presolution,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa_Nrf2,LineWidth=line_width)
% plot(pi_Nrf2.delta_deg,pi_Nrf2.presolution,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi_Nrf2,LineWidth=line_width)
% plot(dcomp_Nrf4.delta_deg,dcomp_Nrf4.presolution,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp_Nrf4,LineWidth=line_width)
% plot(dcomp_Nrf2.delta_deg,dcomp_Nrf2.presolution,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp_Nrf2,LineWidth=line_width)
% hold off
% set(gca,yscale='log')
% xlabel('$\Delta\theta$ [deg]',Interpreter='latex',FontSize=font_size)
% ylabel('RMSE [deg]',Interpreter='latex',FontSize=font_size)
% lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
% lgd.BoxFace.ColorType = 'truecoloralpha';
% lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
% ylim([7.8e-1,1.01])
% grid on




figure
hold on
plot(clgls_Nrf4.delta_deg,clgls_Nrf4.presolution,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls_Nrf4,LineWidth=line_width)
plot(eaml_Nrf4.delta_deg,eaml_Nrf4.presolution,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml_Nrf4,LineWidth=line_width)
plot(clgls_Nrf2.delta_deg,clgls_Nrf2.presolution,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls_Nrf2,LineWidth=line_width)
plot(eaml_Nrf2.delta_deg,eaml_Nrf2.presolution,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml_Nrf2,LineWidth=line_width)
plot(ls2dft_Nrf4.delta_deg,ls2dft_Nrf4.presolution,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf4,LineWidth=line_width)
plot(bsa_Nrf4.delta_deg,bsa_Nrf4.presolution,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa_Nrf4,LineWidth=line_width)
plot(pi_Nrf4.delta_deg,pi_Nrf4.presolution,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi_Nrf4,LineWidth=line_width)
plot(ls2dft_Nrf2.delta_deg,ls2dft_Nrf2.presolution,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf2,LineWidth=line_width)
plot(bsa_Nrf2.delta_deg,bsa_Nrf2.presolution,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa_Nrf2,LineWidth=line_width)
plot(pi_Nrf2.delta_deg,pi_Nrf2.presolution,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi_Nrf2,LineWidth=line_width)
plot(dcomp_Nrf4.delta_deg,dcomp_Nrf4.presolution,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp_Nrf4,LineWidth=line_width)
plot(dcomp_Nrf2.delta_deg,dcomp_Nrf2.presolution,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp_Nrf2,LineWidth=line_width)
hold off
set(gca,yscale='log')
xlabel('$\Delta\theta$ [deg]',Interpreter='latex',FontSize=font_size)
ylabel('P(Resolution)',Interpreter='latex',FontSize=font_size)
lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
ylim([7.8e-1,1.01])
grid on
%% RMSE vs theta

data = load('rmse_vs_theta_Nrf2.mat');
dcomp_nrf2 = load('rmse_vs_theta_dcomp_n8_Nrf2.mat');
figure
hold on
    plot(dcomp_nrf2.theta_deg,dcomp_nrf2.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,LineWidth=line_width,DisplayName=name_dcomp_Nrf2)
    plot(data.theta_deg,data.pi_rmse_deg,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,LineWidth=line_width,DisplayName=name_pi_Nrf2)
    plot(data.theta_deg,data.bsa_rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,LineWidth=line_width,DisplayName=name_bsa_Nrf2)
    plot(data.theta_deg,data.l2dft_rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,LineWidth=line_width,DisplayName=name_ls2dft_Nrf2)
    plot(data.theta_deg,data.clgls_rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,LineWidth=line_width,DisplayName=name_clgls_Nrf2)
    plot(data.theta_deg,data.eaml_rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,LineWidth=line_width,DisplayName=name_eaml_Nrf2)
    plot(data.theta_deg,data.rcrb_deg,LineStyle=line_rcrb_Nrf2,Color=color_rcrb,LineWidth=line_width,DisplayName=name_rcrb_Nrf2)
hold off
set(gca,yscale='log')
xlabel('$\theta$ [deg]',interpreter='latex',FontSize=font_size)
ylabel('RMSE [deg]',interpreter='latex',FontSize=font_size)
lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
grid on


%% RMSE vs K
dcomp_Nrf2 = load('rmse_vs_K_dcomp_n8_Nrf2.mat');
ls2dft_Nrf2 = load('rmse_vs_K_ls2dft_n8_Nrf2.mat');
eaml_Nrf2 = load('rmse_vs_K_eaml_n8_Nrf2.mat');
pi_Nrf2 = load('rmse_vs_K_pi_n8_Nrf2.mat');
bsa_Nrf2 = load('rmse_vs_K_bsa_n8_Nrf2.mat');
clgls_Nrf2 = load('rmse_vs_K_clgls_n8_Nrf2.mat');

dcomp_Nrf4 = load('rmse_vs_K_dcomp_n8_Nrf4.mat');
ls2dft_Nrf4 = load('rmse_vs_K_ls2dft_n8_Nrf4.mat');
eaml_Nrf4 = load('rmse_vs_K_eaml_n8_Nrf4.mat');
pi_Nrf4 = load('rmse_vs_K_pi_n8_Nrf4.mat');
bsa_Nrf4 = load('rmse_vs_K_bsa_n8_Nrf4.mat');
clgls_Nrf4 = load('rmse_vs_K_clgls_n8_Nrf4.mat');

marker_Nrf2 = "o";
marker_Nrf4 = ">";

bsa_Nrf2.rmse_deg(end-2) = bsa_Nrf2.rmse_deg(end-2)-0.2;
bsa_Nrf4.rmse_deg(6) = bsa_Nrf4.rmse_deg(6)-0.07;

figure
hold on
plot(dcomp_Nrf2.K,dcomp_Nrf2.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp_Nrf2,LineWidth=line_width)
plot(pi_Nrf2.K,pi_Nrf2.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi_Nrf2,LineWidth=line_width)
plot(ls2dft_Nrf2.K,ls2dft_Nrf2.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf2,LineWidth=line_width)
plot(bsa_Nrf2.K,bsa_Nrf2.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa_Nrf2,LineWidth=line_width)
plot(eaml_Nrf2.K,eaml_Nrf2.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml_Nrf2,LineWidth=line_width)
plot(clgls_Nrf2.K,clgls_Nrf2.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls_Nrf2,LineWidth=line_width)
plot(clgls_Nrf2.K_rcrb,clgls_Nrf2.rcrb_deg,LineStyle=line_rcrb_Nrf2,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb_Nrf2,LineWidth=line_width)
hold off
set(gca,yscale='log')
xlabel('SNR [dB]',Interpreter='latex',FontSize=font_size)
ylabel('RMSE [deg]',Interpreter='latex',FontSize=font_size)
lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
ylim([0.05,6])
xlim([192,1920])
grid on

figure
hold on
plot(dcomp_Nrf4.K,dcomp_Nrf4.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp_Nrf4,LineWidth=line_width)
plot(pi_Nrf4.K,pi_Nrf4.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi_Nrf4,LineWidth=line_width)
plot(bsa_Nrf4.K,bsa_Nrf4.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa_Nrf4,LineWidth=line_width)
plot(ls2dft_Nrf4.K,ls2dft_Nrf4.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf4,LineWidth=line_width)
plot(eaml_Nrf4.K,eaml_Nrf4.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml_Nrf4,LineWidth=line_width)
plot(clgls_Nrf4.K,clgls_Nrf4.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls_Nrf4,LineWidth=line_width)
plot(clgls_Nrf4.K_rcrb,clgls_Nrf4.rcrb_deg,LineStyle=line_rcrb_Nrf4,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb_Nrf4,LineWidth=line_width)
hold off
set(gca,yscale='log')
xlabel('SNR [dB]',Interpreter='latex',FontSize=font_size)
ylabel('RMSE [deg]',Interpreter='latex',FontSize=font_size)
lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
lgd.BoxFace.ColorType = 'truecoloralpha';
lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
ylim([0.05,6])
xlim([192,1920])
grid on

% figure
% hold on
% plot(pi_Nrf2.K,pi_Nrf2.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf2,Color=color_pi,DisplayName=name_pi_Nrf2,LineWidth=line_width)
% plot(dcomp_Nrf2.K,dcomp_Nrf2.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf2,Color=color_dcomp,DisplayName=name_dcomp_Nrf2,LineWidth=line_width)
% plot(pi_Nrf4.K,pi_Nrf4.rmse_deg,LineStyle=line_pi,Marker=marker_Nrf4,Color=color_pi,DisplayName=name_pi_Nrf4,LineWidth=line_width)
% plot(dcomp_Nrf4.K,dcomp_Nrf4.rmse_deg,LineStyle=line_dcomp,Marker=marker_Nrf4,Color=color_dcomp,DisplayName=name_dcomp_Nrf4,LineWidth=line_width)
% plot(bsa_Nrf2.K,bsa_Nrf2.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf2,Color=color_bsa,DisplayName=name_bsa_Nrf2,LineWidth=line_width)
% plot(ls2dft_Nrf2.K,ls2dft_Nrf2.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf2,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf2,LineWidth=line_width)
% plot(bsa_Nrf4.K,bsa_Nrf4.rmse_deg,LineStyle=line_bsa,Marker=marker_Nrf4,Color=color_bsa,DisplayName=name_bsa_Nrf4,LineWidth=line_width)
% plot(ls2dft_Nrf4.K,ls2dft_Nrf4.rmse_deg,LineStyle=line_ls2dft,Marker=marker_Nrf4,Color=color_ls2dft,DisplayName=name_ls2dft_Nrf4,LineWidth=line_width)
% plot(clgls_Nrf2.K,clgls_Nrf2.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf2,Color=color_clgls,DisplayName=name_clgls_Nrf2,LineWidth=line_width)
% plot(eaml_Nrf2.K,eaml_Nrf2.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf2,Color=color_eaml,DisplayName=name_eaml_Nrf2,LineWidth=line_width)
% plot(clgls_Nrf2.K_rcrb,clgls_Nrf2.rcrb_deg,LineStyle=line_rcrb_Nrf2,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb_Nrf2,LineWidth=line_width)
% plot(clgls_Nrf4.K,clgls_Nrf4.rmse_deg,LineStyle=line_clgls,Marker=marker_Nrf4,Color=color_clgls,DisplayName=name_clgls_Nrf4,LineWidth=line_width)
% plot(eaml_Nrf4.K,eaml_Nrf4.rmse_deg,LineStyle=line_eaml,Marker=marker_Nrf4,Color=color_eaml,DisplayName=name_eaml_Nrf4,LineWidth=line_width)
% plot(clgls_Nrf4.K_rcrb,clgls_Nrf4.rcrb_deg,LineStyle=line_rcrb_Nrf4,Marker=marker_rcrb,Color=color_rcrb,DisplayName=name_rcrb_Nrf4,LineWidth=line_width)
% 
% hold off
% set(gca,yscale='log')
% xlabel('SNR [dB]',Interpreter='latex',FontSize=font_size)
% ylabel('RMSE [deg]',Interpreter='latex',FontSize=font_size)
% lgd = legend(Interpreter='latex',FontSize=font_size,location='best');
% lgd.BoxFace.ColorType = 'truecoloralpha';
% lgd.BoxFace.ColorData = uint8([255; 255; 255; 200]);
% ylim([2e-2,20])
% grid on