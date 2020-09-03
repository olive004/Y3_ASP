% CW2.5

%% 2.5.1

xECG = load('RAW.mat');     % iAMP data
fsECG = 1000;               % read from iAmp_import_v40.m '[fs, res] = read_Header(head_mode)'

channel = 1;                % for some reason our data has 3 channels :3

trial1_cut = 250*fsECG;     % ie Trial 1 ends at 250s
trial2_cut = 500*fsECG;
trial3_cut = 750*fsECG;


xECG_trial1 = xECG.data(1:trial1_cut,channel);
xECG_trial2 = xECG.data(trial1_cut+1:trial2_cut,channel);
xECG_trial3 = xECG.data(trial2_cut+1:trial3_cut,channel);


% Keep as separate vectors cuz not sure if RRI returns same length 
ampthresh1 = 0.1; ampthresh2 = 0.05; ampthresh3 = 0.02;
[xRRI_trial1,fsRRI_trial1] = ECG_to_RRI(xECG_trial1,fsECG);     % removed all anomalies
[xRRI_trial2,fsRRI_trial2] = ECG_to_RRI(xECG_trial2,fsECG);
[xRRI_trial3,fsRRI_trial3] = ECG_to_RRI(xECG_trial3,fsECG);

% load('xRRI_trial1'); load('xRRI_trial2'); load('xRRI_trial3');

% Heartrates
h1 = (ones(size(xRRI_trial1))*60) ./ xRRI_trial1;
h2 = (ones(size(xRRI_trial2))*60) ./ xRRI_trial2;
h3 = (ones(size(xRRI_trial3))*60) ./ xRRI_trial3;

% Get smoothed heartrate with alpha=1
step_size = 10;

alpha = 1;
h1_smooth_1 = smooth_ma(h1, step_size, alpha); 
h2_smooth_1 = smooth_ma(h2, step_size, alpha); 
h3_smooth_1 = smooth_ma(h3, step_size, alpha); 

% Get smoothed heartrate with alpha=0.6
alpha = 0.6;
h1_smooth_06 = smooth_ma(h1, step_size, alpha); 
h2_smooth_06 = smooth_ma(h2, step_size, alpha); 
h3_smooth_06 = smooth_ma(h3, step_size, alpha); 


%% 2.5.1 Plot probability density estimate (PDE) of heart rate vs. h_smooth

figure; hold on
histogram(h1, 'Normalization', 'pdf');
histogram(h2, 'Normalization', 'pdf');
histogram(h3, 'Normalization', 'pdf');
legend('trial 1', 'trial2', 'trial 3');
xlabel('time (s)'); ylabel('h[n] signal');
title('PDE Heartrate');
xlim([0, 200]);
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q251_rri', 'epsc');
hold off

figure; hold on
histogram(h1_smooth_1, 'Normalization', 'pdf');
histogram(h2_smooth_1, 'Normalization', 'pdf');
histogram(h3_smooth_1, 'Normalization', 'pdf');
legend('trial 1', 'trial2', 'trial 3');
xlabel('time (s)'); ylabel('h[n] signal');
xlim([0, 200]);
title('PDE Smooth Heartrate, alpha=1');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q251_1', 'epsc');
hold off

figure; hold on
histogram(h1_smooth_06, 'Normalization', 'pdf');
histogram(h2_smooth_06, 'Normalization', 'pdf');
histogram(h3_smooth_06, 'Normalization', 'pdf');
legend('trial 1', 'trial2', 'trial 3');
xlabel('time (s)'); ylabel('h[n] signal');
xlim([0, 200]);
title('PDE Smooth Heartrate, alpha=0.6');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW2_q251_06', 'epsc');
hold off

%% 2.5.3 Autocorr

% zero mean data using detrend and zscore
xRRI_trial1_z = detrend(xRRI_trial1);
xRRI_trial2_z = detrend(xRRI_trial2);
xRRI_trial3_z = detrend(xRRI_trial3);

xRRI_trial1_z = zscore(xRRI_trial1_z);
xRRI_trial2_z = zscore(xRRI_trial2_z);
xRRI_trial3_z = zscore(xRRI_trial3_z);

% autocorr
[auto1, lag1] = xcorr(xRRI_trial1_z);
[auto2, lag2] = xcorr(xRRI_trial2_z);
[auto3, lag3] = xcorr(xRRI_trial3_z);

% Plot autocorrelations
figure('PaperPosition',[0 0 30 7]); ; hold on
plot(lag1, auto1); plot(lag2, auto2); plot(lag3, auto3);
title('RRI Autocorrelation');
xlabel('time lag tau'); ylabel('Correlation strength');
legend('Trial 1', 'Trial 2', 'Trial 3');
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW2_q253', 'epsc');



%% 2.5.4 AR signal model order

% Stem plot of reflection coefficients up to order p=10
% MDL, AIC, AICc plot (get cum error)
p_max = 10;

N1 = length(xRRI_trial1);
N2 = length(xRRI_trial2);
N3 = length(xRRI_trial3);

error_cum1 = zeros(1, p_max);
error_cum2 = zeros(1, p_max);
error_cum3 = zeros(1, p_max);
MDL1 = zeros(1, p_max);
MDL2 = zeros(1, p_max);
MDL3 = zeros(1, p_max);
AIC1 = zeros(1, p_max);
AIC2 = zeros(1, p_max);
AIC3 = zeros(1, p_max);
AICc1 = zeros(1, p_max);     % Corrected Akaike's Info Crit
AICc2 = zeros(1, p_max);
AICc3 = zeros(1, p_max);

for p=1:p_max
    [part_coeffs1, var1, rc1] = aryule(xRRI_trial1_z, p);
    [part_coeffs2, var2, rc2] = aryule(xRRI_trial2_z, p);
    [part_coeffs3, var3, rc3] = aryule(xRRI_trial3_z, p);
    
    % AR coefficient b must have negative since aryule scales it that way
    % or something
    est_x1 = filter(-part_coeffs1, 1, xRRI_trial1_z);
    est_x2 = filter(-part_coeffs2, 1, xRRI_trial2_z);
    est_x3 = filter(-part_coeffs3, 1, xRRI_trial3_z);
    
    % Cumulative error f
    error_cum1(p) = sum((est_x1 - xRRI_trial1_z).^2);
    error_cum2(p) = sum((est_x2 - xRRI_trial2_z).^2);
    error_cum3(p) = sum((est_x3 - xRRI_trial3_z).^2);

    % MDL
    MDL1(p) = log(error_cum1(p)) + ((p*log(N1)) / N1);
    MDL2(p) = log(error_cum2(p)) + ((p*log(N2)) / N2);
    MDL3(p) = log(error_cum3(p)) + ((p*log(N3)) / N3);
    % AIC
    AIC1(p) = log(error_cum1(p)) + ((2*p) / N1);
    AIC2(p) = log(error_cum2(p)) + ((2*p) / N2);
    AIC3(p) = log(error_cum3(p)) + ((2*p) / N3);
    % AICc
    AICc1(p) = AIC1(p) + ((2*p*(p+1))/(N1-p-1));
    AICc2(p) = AIC2(p) + ((2*p*(p+1))/(N2-p-1));
    AICc3(p) = AIC3(p) + ((2*p*(p+1))/(N3-p-1));
end
% normalizing 
MDL1 = MDL1./max(MDL1); 
AIC1 = AIC1./max(AIC1);
AICc1 = AICc1./max(AICc1);
% normalizing 
MDL2 = MDL2./max(MDL2); 
AIC2 = AIC2./max(AIC2);
AICc2 = AICc2./max(AICc2);
% normalizing 
MDL3 = MDL3./max(MDL3); 
AIC3 = AIC3./max(AIC3);
AICc3 = AICc3./max(AICc3);

p = 1:p_max;

figure('PaperPosition',[0 0 30 6]); 
subplot(1,3,1); hold on
plot(p, MDL1'); plot(p, AIC1'); plot(p, AICc1'); % plot(error(N,:)');
title('Model Order: Trial 1');
xlabel('Model order p'); ylabel('Strength');
lgd = legend('MDL', 'AIC', 'AICc'); %, 'Cumulative sq error');
lgd.Location = 'northwest';
set(gca, 'Fontsize', 18);
hold off
subplot(1,3,2); hold on
plot(p, MDL1'); plot(p, AIC1'); plot(p, AICc1'); % plot(error(N,:)');
title('Model Order: Trial 2');
xlabel('Model order p'); ylabel('Strength');
lgd = legend('MDL', 'AIC', 'AICc'); %, 'Cumulative sq error');
lgd.Location = 'northwest';
set(gca, 'Fontsize', 18);
hold off
subplot(1,3,3); hold on
plot(p, MDL1'); plot(p, AIC1'); plot(p, AICc1'); % plot(error(N,:)');
title('Model Order: Trial 3');
xlabel('Model order p'); ylabel('Strength');
lgd = legend('MDL', 'AIC', 'AICc'); %, 'Cumulative sq error');
lgd.Location = 'northwest';
set(gca, 'Fontsize', 18);
hold off
saveas(gcf,'CW2_q254','epsc');

% 2.5.4 PCF Reflection coeffs
figure('PaperPosition',[0 0 20 8]); hold on
stem(-rc1); stem(-rc2); stem(-rc3);
title('PCF Reflection coefficients iAMP');
xlabel('Model order p'); ylabel('Coefficient');
lgd = legend('Trial 1', 'Trial 2', 'Trial 3'); % lgd.Location = 'southeast';
set(gca, 'Fontsize', 20);
hold off
saveas(gcf,'CW2_q254_pcf','epsc');








%% Checcing w Mings results

% ahahaha sike






