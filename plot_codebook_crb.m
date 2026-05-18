clear,


n = 4;
nrf = 2;
lambda = 1;
d = lambda/2;

theta_deg = linspace(-85,85,1024);
theta = theta_deg*pi/180;

M = ceil(n/(nrf-1));

u_vec = (0:n-1);
F = 1/sqrt(n)*exp(1j*2*pi/n*u_vec.*u_vec.');
codebook = zeros(n,nrf,M);




idx = 0:nrf-1;
for i = 1 : M
    codebook(:,:,i) = F(:,idx+1);
    idx = mod(idx + nrf-1,n);
end

Rs = 10;
sigman2 = 1;


K = 4*10*M;
Km = K/M; 

M = ceil(n/(nrf-1));

rcrb_deg_per_batch = zeros(M,length(theta_deg));
rcrb_deg = zeros(1,length(theta_deg));
% rcrb1_deg = zeros(1,length(theta_deg));


red = [255,0,0]/255;
orange = [255,128,0]/255;
green = [0,255,0]/255;
blue = [0,0,255]/255;

line_width = 2;
font_size = 14;
for i = 1 : length(theta_deg)
    for m = 1 : M
        [CRB,~] = ToolsObj.crb_unc_uncorrelated_theta(d,lambda,theta(i),Rs,sigman2,codebook(:,:,m),Km);
        rcrb_deg_per_batch(m,i) = sqrt(CRB(1,1))*180/pi;
    end
    [CRB,~] = ToolsObj.crb_unc_uncorrelated_theta(d,lambda,theta(i),Rs,sigman2,codebook,Km);
    rcrb_deg(i) = sqrt(CRB(1,1))*180/pi;
end

figure
hold on
plot(theta_deg,rcrb_deg_per_batch(1,:),Color=red,LineWidth=line_width, DisplayName="RCRB only with $\textbf{B}_0$")
plot(theta_deg,rcrb_deg_per_batch(2,:),Color=orange,LineWidth=line_width, DisplayName="RCRB only with $\textbf{B}_1$")
plot(theta_deg,rcrb_deg_per_batch(3,:),Color=green,LineWidth=line_width, DisplayName="RCRB only with $\textbf{B}_2$")
plot(theta_deg,rcrb_deg_per_batch(4,:),Color=blue,LineWidth=line_width, DisplayName="RCRB only with $\textbf{B}_3$")
plot(theta_deg,rcrb_deg,'-k',LineWidth=line_width,DisplayName="RCRB with $\{\textbf{B}_0,\textbf{B}_1,\textbf{B}_2,\textbf{B}_3\}$")
hold off
set(gca,yscale='log')
xlabel('$\theta$ $[\mathrm{deg}]$',Interpreter='latex',FontSize=font_size)
ylabel('$\mathrm{RCRB}$  $\mathrm{[deg]}$',Interpreter='latex',FontSize=font_size)
xlim([-85,85])
ylim([power(10,-1.5),10])
legend(Interpreter="latex",FontSize=14)
grid on

% line_width = 2;
% font_size = 14;
% figure
% hold on
% semilogy(theta_deg,rcrb0_deg,'-b',DisplayName='$\mathrm{Codebook \ without  \ } \textbf{B}_{3}$',LineWidth=line_width)
% semilogy(theta_deg,rcrb1_deg,'-r', DisplayName='$\mathrm{Codebook \ with \ } \textbf{B}_{4}$',LineWidth=line_width)
% hold off
% set(gca,'YScale','log')
% grid on
% xlabel('$\theta$ $[\mathrm{deg}]$',Interpreter='latex',FontSize=font_size)
% ylabel('$\mathrm{RCRB}$  $\mathrm{[deg]}$',Interpreter='latex',FontSize=font_size)
% xlim([-85,85])
% legend(Interpreter="latex",FontSize=14)

