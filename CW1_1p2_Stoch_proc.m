% Q1.2 AND Q1.3 ASP Coursework 
clear all; close all;


% 1.2.1 
M = 100;
N = 100;

x1 = rp1(M, N);
x2 = rp2(M, N);
x3 = rp3(M, N);

% Mean and std of x123
m1 = mean(x1);
m2 = mean(x2);
m3 = mean(x3);

std1 = std(x1);
std2 = std(x2);
std3 = std(x3);

% Plotting vs. t
figure; hold on
plot(m1); plot(m2); plot(m3); 
title(['Mean of rp, M=', num2str(M)]); xlabel('t'); ylabel('Mean val');
legend('Mean 1', 'Mean 2', 'Mean 3');
set(gca, 'Fontsize', 20);
hold off
plot_name = strcat('CW1_q121_mean_', num2str(M));
saveas(gcf,plot_name,'epsc');

figure; hold on
plot(std1); plot(std2); plot(std3); 
title(['Std of rp, M=', num2str(M)]); xlabel('t'); ylabel('Std val');
legend('Std 1', 'Std 2', 'Std 3'); 
set(gca, 'Fontsize', 20);
hold off
plot_name = strcat('CW1_q121_std_', num2str(M));
saveas(gcf,plot_name,'epsc');



% 1.2.2
M = 1000;   % should be 4, made it big for 1.2.3
N = 1000;

x1 = rp1(M, N);
x2 = rp2(M, N);
x3 = rp3(M, N);

% Mean and std of x123
% m1 = mean(x1, 2);
% m2 = mean(x2, 2);
% m3 = mean(x3, 2);
m1 = mean(mean(x1, 2));
m2 = mean(mean(x2, 2));
m3 = mean(mean(x3, 2));
% std's calc'd with transposed x signal cuz std() acts on first non 1 sized axis
% std1 = std(x1.');     
% std2 = std(x2.');
% std3 = std(x3.');
std1 = mean(std(x1.'));     
std2 = mean(std(x2.'));
std3 = mean(std(x3.'));

% Plotting M=4 vs. t
figure; hold on
plot(m1); plot(m2); plot(m3); 
title(['Mean of rp, M=', num2str(M)]); xlabel('t'); ylabel('Mean val');
legend('Mean rp1', 'Mean rp2', 'Mean rp3')
hold off

figure; hold on
plot(std1); plot(std2); plot(std3); 
title(['Std of rp, M=', num2str(M)]); xlabel('t'); ylabel('Std val');
legend('Std rp1', 'Std rp2', 'Std rp3')
hold off





% 1.2.3 Mathematical description of each of the three stochastic processes. 
% Calculate the theoretical mean and variance for each case and compare the 
% theoretical results with those obtained by sample averaging.
% --> See the 3 functions for my derivation of the theoretical values

% m_theo1 = 0.02;
% m_theo2 = 7/12;
% m_theo3 = 0.5;
% 
% std_theo1 = 0.25/3 * b^2 (sin( pi/ N))^2 - 0.25;
% std_theo2 = 0.19305;
% std_theo3 = 0.5;







