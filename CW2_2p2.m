% ASP CW Q 2.2
clear all;
close all; 

% q221
N = 1000;
a = [1];
filt_order = 9;
b = ones(filt_order,1);

x = randn(1, N);
y=filter(b,a,x);

[acf, tau] = xcorr(x,y, 'unbiased');

stem(tau, acf);
title(['\fontsize{16}CCF of x and y stem plot, N=', num2str(N)]); 
xlabel('tau'); ylabel('Correlation strength')
xlim([-20 20]);
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q221', 'epsc');