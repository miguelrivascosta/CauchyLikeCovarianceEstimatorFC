clear, close all

N = 32;
Nrf = 8;
K = 10*ceil(N/(Nrf-1));

L_max = 10;
snrdb = 5;

eaml = SheinvaldMethod(N,Nrf);
clgls = CLGLS(N,Nrf);
dcomp = Dcomp(N,Nrf,K,1000);

time_eaml = zeros(1,L_max);
time_clgls = zeros(1,L_max);
time_dcomp = zeros(1,L_max);

for l = 1 : L_max
    theta = deg2rad(linspace(-60,60,l));
    A = exp(1j*pi*(0:N-1)'.*sin(theta));
    Rs = eye(length(theta));
    x = ToolsObj.generate_uncorrelated_gaussian_channel(A,Rs,10^(-snrdb/10),K);

    time_eaml(l) = timeit(@() eaml.estimate_doas(x,l));
    time_clgls(l) = timeit(@() clgls.estimate_doas(x,l));
    time_dcomp(l) = timeit(@() dcomp.estimate_doas(x,l));
    
    fprintf("Iteration: %d/%d\n",l,L_max)
end

figure
hold on
plot(1:L_max,time_eaml,'-o')
plot(1:L_max,time_clgls,'-s')
plot(1:L_max,time_dcomp,'->')
hold off
set(gca,yscale="log")