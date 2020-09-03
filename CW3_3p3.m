% ASP CW 3.3
clear all
close all

%% 3.3.3 Sunspot LSE: get AR coeffs

load('sunspot.dat')
time = sunspot(:,1);
sunspot = sunspot(:,2);
sun_norm = zscore(sunspot);
N = length(sun_norm);

rxx = xcorr(sun_norm,  'biased'); %, N); %
offset = 0; % 2; %
rxx_cent = N+offset/2;
in_sample = rxx(rxx_cent:end-1);
out_sample = rxx(rxx_cent+1:end);

x = in_sample;

% a = (H^T H)^-1 H x
p_max = 10;
a = zeros(p_max,p_max);

% Get MDL, AIC, & AICc
J = zeros(1,p_max);
MDL = zeros(1, p_max);      
AIC = zeros(1, p_max);
AICc = zeros(1, p_max);     % Corrected Akaike's Info Crit
M = N-1;  
for p=1:p_max
    % Get H
    r = in_sample(1:M);
    c = in_sample(1:p);
    H = toeplitz(r,c); % toeplitz(rxx(N+offset/2:end));
%     H = H(1:M, 1:p);
%     H = zeros(M, p);
%     for k = 1:M
%         for i=1:p
%             H(k, i) = rxx(rxx_cent+k - i);
%         end
%     end
    
    % Get a (aka theta_hat)
    theta_hat = ((H' * H) \ H') * out_sample;
    a(p,:) = [theta_hat', zeros(1, p_max - p)];
    
    % J_min: Get error by plugging theta_hat into J
    s = H * a(p, 1:p)';
    error = (out_sample - s); 
    J(p) = (error'*error);                      % Minimum J occurs at theta_hat
    
    MDL(p) = log(J(p)) + (p*log(N))/N;
    AIC(p) = log(J(p)) + (2*p)/N;
    AICc(p) = AIC(p) + ((2*p*(p+1))/(N-p-1));
end
% normalizing 
MDL = MDL./max(MDL); 
AIC = AIC./max(AIC);
AICc = AICc./max(AICc);

p = 3;
% sys = ar(sun_norm, p, 'ls');
% est_x = filter(1, sys.A(1:p), sun_norm);
est_x = filter(1,a(p,1:p), sun_norm); % sunspot);
est_x = zscore(est_x);

figure('PaperPosition', [0 0 30 10]); subplot(1,2,1); hold on  
plot(time, est_x); 
plot(time, sun_norm); % sunspot); 
title(['LSE AR model p=', num2str(p)]);
xlabel('time (years)'); ylabel('Sunspots');
set(gca, 'Fontsize', 18);
legend('AR(p)','sunspot signal'); 
hold off

subplot(1,2,2); hold on;
p = 1:p_max;
plot(p, MDL'); plot(p, AIC'); plot(p, AICc'); % plot(error(N,:)');
title('Model Order Estimation');
xlabel('Model order p'); ylabel('Strength');
set(gca, 'Fontsize', 18);
legend('MDL', 'AIC', 'AICc'); %, 'Cumulative sq error');
hold off
saveas(gcf,'CW3_q333','epsc');


%% 3.3.5 Power spectra

p=3;

sys = ar(sun_norm, p, 'ls');
a_p = [a(p,1:p)];       % negative cuz of freqz definition
% a_p = sys.A(1:p);
psd_sun = pgm(sun_norm);
psd_sun = psd_sun ./ max(psd_sun);
[h,w]=freqz(1,a_p,N);
est_sun = (abs(h).^2) / max(abs(h).^2);
f = 0:(1/N):((N/2-1)/N);

figure('PaperPosition', [0 0 25 10]); subplot(1,2,1);
hold on
plot(f, 10*log10(psd_sun(1:N/2))); 
plot(w/(2*pi), 10*log10(est_sun),'b', 'LineWidth', 2);       % weird 8*abs(h).^2/N
title(['PSD Comparison, p=', num2str(p)]);
xlabel('Normalized frequency (Hz)'); ylabel('Power (dB)');
lgd = legend('Original PSD', 'Estimated PSD'); lgd.Location = 'southeast';
set(gca, 'Fontsize', 20);
hold off
subplot(1,2,2);
hold on
plot(f, (psd_sun(1:N/2))); 
plot(w/(2*pi), (est_sun),'b', 'LineWidth', 2);       % weird 8*abs(h).^2/N
title(['PSD non-dB, p=', num2str(p)]);
xlabel('Normalized frequency (Hz)'); ylabel('Power');
lgd = legend('Original PSD', 'Estimated PSD'); lgd.Location = 'northeast';
set(gca, 'Fontsize', 20);
hold off

saveas(gcf, 'CW3_q335', 'epsc');




%% 3.3.6 MSE of J

rxx = xcorr(sun_norm, 'biased'); %  N); ,
offset = 0;
in_sample = rxx(rxx_cent:end-1);
out_sample = rxx(rxx_cent+1:end);

x = in_sample;

N = [10 : 5 : 250];
NN = length(N);
N_s = length(sunspot);

p = 3;

% Get MDL, AIC, & AICc
error_cum = zeros(1, NN);
J = zeros(1,NN);

i = 1;
for n = N
    
    % Get H
    r = in_sample(1:n);
    c = in_sample(1:p);
    H = toeplitz(r,c); % toeplitz(rxx(N+offset/2:end));
%     H = H(1:M, 1:p);
%     H = zeros(M, p);
%     for k = 1:M
%         for i=1:p
%             H(k, i) = rxx(rxx_cent+k - i);
%         end
%     end
    
    % Get a (aka theta_hat)
    theta_hat = ((H' * H) \ H') * out_sample(1:n);
    a(p,:) = [theta_hat', zeros(1, p_max - p)];
    
    % J_min: Get error by plugging theta_hat into J
    s = H * a(p, 1:p)';
    error = (out_sample - s); 
    J(p) = (error'*error);                      % Minimum J occurs at theta_hat
    
    i = i+1;
end
% normalizing 
MSE = MSE./max(MSE);

figure; hold on
plot(N, MSE); %plot(N, J);
title(['MSE']);
xlabel('Data length N'); ylabel('MSE');
legend('MSE', 'Normalized J');
set(gca, 'Fontsize', 24);
hold off
saveas(gcf, 'CW3_q336', 'epsc');

