close all

%% 
clear

blue="#0072BD";
orange="#D95319";
yellow="#EDB120";

data0 = load("data_complexity_trace0.mat");
data1 = load("data_complexity_trace1.mat");
data2 = load("data_complexity_woodbury_inversion_var_nrf.mat");

lw_measured = 1.5;
lw_fit = 2.5;
font_size_legend =  14;
font_size_xlabel = 14;
font_size_ylabel = 14;
figure
hold on

plot(data2.log2_nrf,data2.log2_T,'o',Color=yellow,LineWidth=lw_measured,DisplayName='$ T(\mathbf{\hat{\mathbf{s}}_{\mathrm{cl,e}}})$')
plot(data2.log2_nrf_interp,data2.log2_T_interp,'-k',LineWidth=lw_fit,HandleVisibility="off")

plot(data1.log2_nrf,data1.log2_T,'o',Color=blue,LineWidth=lw_measured,DisplayName='$T(\mathbf{\Sigma}_{\mathrm{cl}}^H\hat{\mathbf{R}}_{\epsilon,m}^{-1}\hat{\mathbf{p}})$')
plot(data1.log2_nrf_interp,data1.log2_T_interp,'-k',LineWidth=lw_fit,HandleVisibility="off")

plot(data0.log2_nrf,data0.log2_T,'o',Color=orange,LineWidth=lw_measured,DisplayName='$ T(\mathbf{\Sigma}_{\mathrm{cl}}^H\hat{\mathbf{R}}_{\epsilon,m}^{-1}\mathbf{\Sigma}_{\mathrm{cl}})$')
plot(data0.log2_nrf_interp,data0.log2_T_interp,'-k',LineWidth=lw_fit,DisplayName='Power-Law fit')

hold off
xlabel('$\mathrm{log}_{2}(N_{\mathrm{RF}})$','Interpreter','latex',FontSize=font_size_xlabel)
ylabel('$\mathrm{log}_2(\mathrm{Time})$','Interpreter','latex',FontSize=font_size_ylabel)
grid on

lgd = legend(FontSize=font_size_legend);
set(lgd, 'Interpreter', 'latex');
ylim([-6 4])
xlim([7.5 12.8])
