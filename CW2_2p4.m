% Quetion 2.4 Omg
clear all
close all
clc

%% q2.4.1
load('NASDAQ.mat');
closing = NASDAQ.Close;
date = NASDAQ.Date;
N = length(closing);

p_max = 10;         % Largest order to analyze, we only care about AR(1)
closing_z = zscore(closing);
[ar1_coeff, ar1_var, ar1_rc] = aryule(closing_z, p_max);
ar1 = filter(-ar1_coeff, 1, closing);

% Plot of closing prices and AR(1) closing
% figure; hold on
% plot(date, closing); plot(date, ar1); 
% title('NASDAQ AR(1) Closing');
% xlabel('Date'); ylabel('Price');
% lgd = legend('Og. closing', 'AR(1) closing');
% lgd.Location = 'southeast';
% set(gca, 'Fontsize', 20);
% saveas(gcf, 'CW2_q241', 'epsc');
% hold off

% Plot of Partial Correlation
figure('PaperPosition',[0 0 30 8]);
subplot(1,2,1); hold on
stem(-ar1_rc); plot(zeros(1,p_max)+ar1_var); plot(zeros(1,p_max)-ar1_var);
lgd = legend('Ref coeffs', 'upper bound', 'lower bound');
set(gca, 'Fontsize', 20);
title('AR Reflection Coefficients'); xlabel('Model order'); ylabel('RC')
hold off

% Plot MDL and AIC
error_cum = zeros(1, p_max);
MDL = zeros(1, p_max);      
AIC = zeros(1, p_max);

for p = 1:p_max
    error_cum(p) = sum((ar1 - closing_z).^2);
    MDL(p) = log(error_cum(p)) + ((p*log(N)) / N);
    AIC(p) = log(error_cum(p)) + ((2*p) / N);
end
% normalizing 
MDL = MDL./max(MDL); 
AIC = AIC./max(AIC);

subplot(1,2,2); hold on
plot(MDL); plot(AIC);
title('NASDAQ MDL and AIC'); 
xlabel('Model Order'); ylabel('Measure Strength');
lgd = legend('MDL', 'AIC'); lgd.Location = 'northwest';
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW2_q241_AIC_rc', 'epsc');



%% 2.4.3
N = [1 : 50 : 1001];
true_sig_sq = [1 : 50 : 1001];
[N_m, true_sig_sq_m] = meshgrid(N, true_sig_sq);    % turns each vector into matrix
rxx = xcorr(closing_z, 'unbiased');

% Get CRLB of I(theta) parameters
crlb_a = true_sig_sq_m ./ (N_m.*(rxx(length(closing))));
crlb_sig_sq = (2*true_sig_sq_m.^2) ./ N_m;

% Heatmap of CRLB 
figure('PaperPosition',[0 0 30 8]); subplot(1,2,1);
h = heatmap(N, true_sig_sq, crlb_a);
h.Colormap = jet;
h.ColorScaling = 'log';
title('\fontsize{20}Estimatetd sigma^2 CRLB Heatmap');
xlabel('\fontsize{20}N'); ylabel('\fontsize{20}Variance');

subplot(1,2,2);
h = heatmap(N, true_sig_sq, crlb_sig_sq);
h.Colormap = jet;
h.ColorScaling = 'log';
title('\fontsize{20}Estimatetd a1 CRLB Heatmap');
xlabel('\fontsize{20}N'); ylabel('\fontsize{20}Variance');

saveas(gcf, 'CW2_q243', 'epsc');


%% 2.4.1cii

N = length(ar1);
crlb_a1_ar1 = (1/N) * (1 - (ar1_rc(1)).^2);
a1_unity = linspace(ar1_rc(1), -1, N);
crlb_a1 = (1/N) * (1 - a1_unity.^2);

figure('PaperPosition',[0 0 10 8]); plot(a1_unity, crlb_a1);
xlabel('a_1'); ylabel('CRLB of a_1');
title('CRLB vs. a_1');
text(a1_unity(2*length(a1_unity)/3), crlb_a1(3*length(crlb_a1)/4), ['\fontsize{20} CRLB from AR(1): ', num2str(crlb_a1_ar1)]);
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q243ii', 'epsc');









