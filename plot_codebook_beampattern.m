clear, close all


n = 4;
nrf = 2;

F = 1/sqrt(n)*exp(1j*2*pi/n*(0:n-1).*(0:n-1).');

if n==nrf
	M = 1;
else
	M = ceil(n/(nrf-1));
end

configuration_matrix = zeros(M,nrf);
configuration_matrix(1,:) = 0:nrf-1;
for i = 2 : M
	configuration_matrix(i,:) = mod(configuration_matrix(i-1,:)+nrf-1,n);
end

codebook = zeros(n,nrf,M);
for i = 1 : M
	codebook(:,:,i) = F(:,configuration_matrix(i,:)+1);
end


theta_plot_deg = linspace(-90,90,1000);
theta_plot = theta_plot_deg*pi/180;

a_plot = exp(1j*pi*sin(theta_plot).*(0:n-1).');
full_beam_pattern = abs(F'*a_plot).^2;





logarithmic = 0;
x_limits = [-90,90];
font_size = 14;
line_width = 2;
red = [255,0,0]/255;
orange = [255,128,0]/255;
green = [0,255,0]/255;
blue = [0,0,255]/255;

if logarithmic == 0
	G_full = (abs(F'*a_plot).^2);
	beampattern_B0 = ((sum(abs(codebook(:,:,1)'*a_plot).^2)));
	beampattern_B1 = ((sum(abs(codebook(:,:,2)'*a_plot).^2)));
	beampattern_B2 = ((sum(abs(codebook(:,:,3)'*a_plot).^2)));
	beampattern_B3 = ((sum(abs(codebook(:,:,4)'*a_plot).^2)));
	y_limits = [0,4.5];

else
	G_full = 10*log10(abs(F'*a_plot).^2);
	beampattern_B0 = 10*log10((sum(abs(codebook(:,:,1)'*a_plot).^2)));
	beampattern_B1 = 10*log10((sum(abs(codebook(:,:,2)'*a_plot).^2)));
	beampattern_B2 = 10*log10((sum(abs(codebook(:,:,3)'*a_plot).^2)));
	beampattern_B3 = 10*log10((sum(abs(codebook(:,:,4)'*a_plot).^2)));
	y_limits = [0.4,30];
end

	subplot(4,1,1)
	hold on
		plot(theta_plot_deg,beampattern_B0,'-',Color=red,LineWidth=line_width)
		% plot(theta_plot_deg,G_full(configuration_matrix(1,:)+1,:),'-k')
	hold off
	xlabel('$\theta$ [deg]',Interpreter='latex',FontSize=font_size)
	% ylabel('$||\textbf{B}_0^H \textbf{a}(\theta)||^2_2$',Interpreter='latex',FontSize=font_size)
	grid on
	% set(gca,'YScale','log')
	% xlim(x_limits)
	ylim(y_limits)
	% set(get(gca,'ylabel'),'rotation',0)
	
	subplot(4,1,2)
	hold on
	% plot(theta_plot_deg,G_full,'--k')
		plot(theta_plot_deg,beampattern_B1,'-',Color=orange,LineWidth=line_width)
		% plot(theta_plot_deg,G_full(configuration_matrix(2,:)+1,:),'-k')
	hold off
	xlabel('$\theta$ [deg]',Interpreter='latex',FontSize=font_size)
	% ylabel('$||\textbf{B}_1^H \textbf{a}(\theta)||^2_2$',Interpreter='latex',FontSize=font_size)
	grid on
	% set(gca,'YScale','log')
	% xlim(x_limits)
	ylim(y_limits)
	
	subplot(4,1,3)
	hold on
		plot(theta_plot_deg,beampattern_B2,'-',Color=green,LineWidth=line_width)
		% plot(theta_plot_deg,G_full(configuration_matrix(3,:)+1,:),'-k')
	hold off
	xlabel('$\theta$ [deg]',Interpreter='latex',FontSize=font_size)
	% ylabel('$||\textbf{B}_2^H \textbf{a}(\theta)||^2_2$',Interpreter='latex',FontSize=font_size)
	grid on
	% set(gca,'YScale','log')
	% set(get(gca,'ylabel'),'rotation',0)
	% xlim(x_limits)
	ylim(y_limits)
	
	subplot(4,1,4)
	hold on
		plot(theta_plot_deg,beampattern_B3,'-',Color=blue,LineWidth=line_width)
		% plot(theta_plot_deg,G_full(configuration_matrix(4,:)+1,:),'-k')
	hold off
	xlabel('$\theta$ [deg]',Interpreter='latex',FontSize=font_size)
	% ylabel('$||\textbf{B}_3^H \textbf{a}(\theta)||^2_2$',Interpreter='latex',FontSize=font_size)
	grid on
	% set(gca,'YScale','log')
	% set(get(gca,'ylabel'),'rotation',0)
	% xlim(x_limits)
	ylim(y_limits)



% for i = 1 : M 
% 	figure
% 	Pm = abs(codebook(:,:,i)'*a_plot).^2;
% 	hold on
% 	plot(theta_plot,10*log10(sum(Pm)),LineWidth=1.5)
% 	plot(theta_plot,10*log10(Pm),'--k',LineWidth=1.5)
% 	hold off
% 	% plot(theta_plot,20*log10(Pm),'--k')
% end
% hold off
% ylim([-45,20])




