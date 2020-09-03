% ASP CW q 2
% q211
x = randn(1, 1000);
[acf, tau] = xcorr(x,'unbiased');

figure('PaperPosition', [0 0 25 10]); hold on
subplot(1,2,2); plot(tau, acf); 
xlim([-999 999]);
title(['Autocorrelation of x']); 
xlabel('tau'); ylabel('Correlation strength');
set(gca, 'Fontsize', 20);
subplot(1,2,1); plot(x);
title(['WGN']); 
xlabel('n'); ylabel('Signal strength');
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW2_q211', 'epsc');


%% q212
figure('PaperPosition', [0 0 10 8]); plot(tau, acf); 
xlim([-50 50]);
title(['Autocorrelation of x zoomed in']); 
xlabel('tau'); ylabel('Correlation strength');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q212', 'epsc');


%% q214
N = 1000;
x = randn(1, N);
a = [1];
filt_order = 9;
b = ones(filt_order,1);
y=filter(b,a,x);
[acf_y, tau_y] = xcorr(y, 'unbiased');

figure('PaperPosition', [0 0 15 12]);
stem(tau_y, acf_y);
title(['MA-filtered autocorrelation stem plot, N=', num2str(N)]); 
xlabel('tau'); ylabel('Correlation strength')
xlim([-20 20]);
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q214', 'epsc');


