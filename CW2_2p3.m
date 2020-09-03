% CW2 q 2.3
clear all;
close all;

% q231
M = 1;
N = 1000;
a1_range = [-2.5, 2.5];
a2_range = [-1.5, 1.5];

a0 = 1;
a1 = a1_range(1) + (a1_range(2) - a1_range(1)) * rand(M, N);
a2 = a2_range(1) + (a2_range(2) - a2_range(1)) * rand(M, N);

N2 = 1000;
w = randn(M, N2);
b = [1];
a = [a0, a1, a2];
x=filter(b,a,w);

figure; hold on
plot(a1, a2, 'r*');
cond1 = (a1 < (1 - a2));
cond2 = (a1 > (a2 -1));
cond3a = a2<1;
cond3b = a2>(-1);
a1_bound = a1( cond1 & cond2 & cond3a & cond3b);
a2_bound = a2( cond1 & cond2 & cond3a & cond3b);
% a2(a1 < (1 - a2));
% a2_bound = a2_bound(a1 > (a2 -1));
% a2_bound = a2_bound( a2<1 );
% a2_bound = a2_bound( a2>(-1));
plot(a1_bound, a2_bound, 'g*');
title(['\fontsize{18}Region of convergence, N=', num2str(N)]);
xlabel('\fontsize{18}a1'); ylabel('\fontsize{18}a2');
legend('unbounded', 'converging')
hold off


%% q232 Sunspot
% Use load sunspot.dat -- The data are contained in the second column of 
% the variable sunspot
load sunspot.dat        % 'sunspot'
sun = sunspot(:,2);
sun_norm = zscore(sun);
time = sunspot(:,1);


% ACF sunspots
figure('PaperPosition', [0 0 35 14]);
plotcount=1;
m_sun = zeros(1,3);
for N = [5, 20, 250]
    sample_sun = sun(1:N);
    [acf, tau] = xcorr(sample_sun, 'unbiased');    
    sample_sun = sun(1:N);
    [acf_m, tau_m] = xcorr(sample_sun - mean(sample_sun), 'unbiased');
    subplot(1,3,plotcount);
    stem(tau, acf); hold on
    stem(tau_m, acf_m);
    title(['\fontsize{18}ACF of sunspots, N=', num2str(N)])
    xlabel('Time lag tau'); ylabel('ACF strength');
    set(gca, 'Fontsize', 18);
    lgd = legend('Sunspot data', 'Zero-mean data'); lgd.Location = 'south';

    plotcount = plotcount+1;
    hold off
end
saveas(gcf,'CW2_q232','epsc');
%%
p =3;
a = aryule(sun_norm, p);
est_x = filter(1, a, sun_norm);
est_x = zscore(est_x);
figure('PaperPosition', [0 0 10 8]); hold on
plot(time, sun_norm); plot(time, est_x);
title('YW AR(2) Sunspot estimate');
xlabel('Time (yrs)'); ylabel('Sunspots');
lgd = legend('Sunspot', 'Estimate'); lgd.FontSize = 15;
set(gca, 'Fontsize', 18);
saveas(gcf, 'CW2_q232_est', 'epsc');




%% q233: YW model order

p_max = 10;
sun_norm = zscore(sunspot(:,2));    % Normalizing sunspot to 0 mean 1 variance
% Calculate AR(p) models 
for p=1:p_max
    [part_coeffs, var, ref_coeff] = aryule(sunspot(:,2), p);
    [part_coeffs_norm, var_norm, ref_coeff_norm] = aryule(sun_norm, p);  
    
end

% Plot stem of 'reflection coefficients'
figure('PaperPosition',[0 0 35 8]); subplot(1,2,1); hold on
stem(-ref_coeff); plot(zeros(1,p)+var_norm); plot(zeros(1,p)-var_norm);
set(gca, 'Fontsize', 20);
title('RCs'); xlabel('Model order');
hold off
subplot(1,2,2); hold on
stem(-ref_coeff_norm); plot(zeros(1,p)+var_norm); plot(zeros(1,p)-var_norm); 
set(gca, 'Fontsize', 20); 
title('RCs normalized'); ylabel('Reflection coeff');
hold off
saveas(gcf,'CW2_q233','epsc');


%% q 2.3.4: Model order estimation: MDL, AIC, AICc
N = length(sunspot);
wgn = randn(N, 1);

p_max = 10;                 % Max. order to estimate to 
error_cum = zeros(1, p_max);
MDL = zeros(1, p_max);      
AIC = zeros(1, p_max);
AICc = zeros(1, p_max);     % Corrected Akaike's Info Crit

% Get MDL, AIC, & AICc
for p=1:p_max
    part_coeffs = aryule(sun_norm, p);
    a = [part_coeffs(2:end)];
    % est_x = filter(-part_coeffs, 1, sun_norm); %filter(1, a, sun_norm);  %
    est_x = filter(1, a, sun_norm);
    
    % Cumulative error f
    error_cum(p) = sum((est_x - sun_norm).^2);

    MDL(p) = log(error_cum(p)) + ((p*log(N)) / N);   
    AIC(p) = log(error_cum(p)) + ((2*p) / N);
    AICc(p) = AIC(p) + ((2*p*(p+1))/(N-p-1));
    
end
% normalizing 
max_i = max([max(MDL), max(AIC), max(AICc)]);
MDL = MDL./max(MDL); 
AIC = AIC./max(AIC);
AICc = AICc./max(AICc);
% MDL = MDL./max_i; 
% AIC = AIC./max_i;
% AICc = AICc./max_i;

a = aryule(sun_norm, 3);
est_x = filter([1], a, sun_norm);
est_x = zscore(est_x);
figure; hold on
plot(est_x); plot(sun_norm);
hold off

figure; hold on;
p = 1:p_max;
plot(p, MDL'); plot(p, AIC'); plot(p, AICc'); % plot(error(N,:)');
title('Model Order Estimation');
xlabel('Model order p'); ylabel('Strength');
set(gca, 'Fontsize', 18);
legend('MDL', 'AIC', 'AICc'); %, 'Cumulative sq error');
hold off
saveas(gcf,'CW2_q234','epsc');



%% 2.3.5 Predicting with AR models

ar1_sys = ar(sun_norm, 1, 'yw');
ar2_sys = ar(sun_norm, 2, 'yw');
ar10_sys = ar(sun_norm, 10, 'yw');


% Forecasting
figure('PaperPosition',[0 0 35 30]);    % [left bottom width height]
m_all = [1,2,5,10];     % prediction horizons
error_ar1 = zeros(1, length(m_all));
error_ar2 = zeros(1, length(m_all));
error_ar10 = zeros(1, length(m_all));
plot_count = 1;
for m = m_all
    est_ar1 = predict(ar1_sys, sun_norm, m);
    est_ar2 = predict(ar2_sys, sun_norm, m);
    est_ar10 = predict(ar10_sys, sun_norm, m);

    error_ar1(plot_count) = sum((est_ar1 - sun_norm).^2);
    error_ar2(plot_count) = sum((est_ar2 - sun_norm).^2);
    error_ar10(plot_count) = sum((est_ar10 - sun_norm).^2);
    
    % plot sunspot + estimate
    subplot(2,2,plot_count); hold on
    plot(sun_norm); 
    plot(est_ar1); plot(est_ar2); plot(est_ar10);
    title(['Forecast, steps m=', num2str(m)]);
    xlabel('time step n'); ylabel('signal strength');
    legend('sunspot', 'AR(1)', 'AR(2)', 'AR(10)');
    set(gca, 'Fontsize', 20);
    hold off
    
    plot_count = plot_count+1;
    % plot error
end   
saveas(gcf,'CW2_q235','epsc');

figure; hold on 
plot(m_all, error_ar1); plot(m_all, error_ar2); plot(m_all, error_ar10); 
title('Error of AR model forecasts');
xlabel('Prediction horizon m'); ylabel('Cumulative squared error');
xticks(m_all);
set(gca, 'Fontsize', 20);
legend('AR(1)', 'AR(2)', 'AR(10)');
hold off
saveas(gcf,'CW2_q235_error','epsc');


