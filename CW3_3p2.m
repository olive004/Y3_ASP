% CW3 3.2
close all; clear all; 

%% 3.2.1 and 3.2.2

N = 1064;
N_cut = N - 40;
b = 1;
a_og = [1, 0.9];       % AR coeffs
cut_i = 40;

x = randn(1,N);
y = filter(b,a_og,x);
y = y((cut_i+1):end);

% PSD
[h,w]=freqz(b,a_og,512);
y_psd = pgm(y);
f = (cut_i/N):(1/N):((N-1)/N);

% Zoom region
f_zoom = [0.4, 0.5];

% Plotting
figure('PaperPosition', [0 0 30 10]); subplot(1,2,1); hold on
plot(f, y_psd, 'r');
plot(w/(2*pi),abs(h).^2,'b', 'LineWidth', 2);
title(['Filtered PSD']);
xlabel('Normalized frequency'); ylabel('Power');
legend('Exact PSD', 'Periodogram');
set(gca,'Xlim',[0,0.5]);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off
subplot(1,2,2); hold on
plot(f, y_psd, 'r'); plot(w/(2*pi),abs(h).^2,'b', 'LineWidth', 2);
title(['Zoomed Filtered PSD']);
xlabel('Normalized frequency'); ylabel('Power');
legend('Exact PSD', 'Periodogram');
set(gca,'Xlim',f_zoom);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off
saveas(gcf, 'CW3_q321', 'epsc');



%% 3.2.4

R_y = xcorr(y, 'unbiased');
a1_hat = -R_y(N_cut + 1) / R_y(N_cut);
sigma2_hat =  R_y(N_cut) + a1_hat* R_y(N_cut + 1);

b = sigma2_hat;
a_og = [1 a1_hat];

% P_y = sigma2 / (abs(1 + a1*exp(- 2*pi*f))^2;
[P_y,w_model] = freqz(b,a_og, N_cut/2); 

% Plotting
figure('PaperPosition', [0 0 30 10]); subplot(1,2,1); hold on
plot(f, y_psd, 'r', 'DisplayName', 'Periodogram');
plot(w/(2*pi),abs(h).^2,'b', 'LineWidth', 2, 'DisplayName', 'Exact PSD');
plot(w_model/(2*pi),abs(P_y).^2,'g', 'LineWidth', 2, 'DisplayName', 'Model PSD');
title(['Model-based PSD']);
xlabel('Normalized frequency (Hz)'); ylabel('Power');
lgd = legend('show'); lgd.Location = 'northwest';
set(gca,'Xlim',[0,0.5]);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off
subplot(1,2,2); hold on
plot(f, y_psd, 'r', 'DisplayName', 'Periodogram');
plot(w/(2*pi),abs(h).^2,'b', 'LineWidth', 2, 'DisplayName', 'Exact PSD');
plot(w_model/(2*pi),abs(P_y).^2,'g', 'LineWidth', 2, 'DisplayName', 'Model PSD');
title(['Zoomed Model-based PSD']);
xlabel('Normalized frequency (Hz)'); ylabel('Power');
lgd = legend('show'); lgd.Location = 'northwest';
set(gca,'Xlim',f_zoom);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off
saveas(gcf, 'CW3_q324', 'epsc');



%% 3.2.5 Sunspot model

load sunspot.dat
% date = sunspot(:,1);
sunspot = sunspot(:, 2);
N = length(sunspot);
f_zoom = [0 0.15];
model_orders = [1,2,5,10];

% Generate PSD's (couldn't be fucked for a for loop)
% original series
sun_og_psd = pgm(sunspot);
Rss = xcorr(sunspot, 'unbiased');       % autocorr og sunspot
% model orders 1,2,5,10
a_og = zeros(1, length(model_orders));
var_og = zeros(1, length(model_orders));
h_og = zeros(N/2, length(model_orders));
w_og = zeros(N/2, length(model_orders));
j = 1;
for i=model_orders
    % Method like above with xcorr
%     a_og(j) = - Rss(N+i) / Rss(N+i-1);
%     var_og(j) = Rss(N+i-1) + a_og(j)*Rss(N+i);
%     
%     a = [1 a_og];
%     b = var_og(j);
%     
%     [h_og(:,j), w_og(:,j)] = freqz(b,a,N/2);
%     j = j+1;
    
    % Aryule method
    [a_og_yule, var_og_yule] = aryule(sunspot, model_orders(j));
    [h_og_yule(:,j), w_og_yule(:,j)] = freqz(var_og_yule^(1/2), a_og_yule, N/2);
    
    j = j+1;
end


% mean-centred (mc) version
sun_mc = zscore(sunspot);
sun_mc_psd = pgm(sun_mc);

Rss = xcorr(sun_mc, 'unbiased');       % autocorr og sunspot
% model orders 1,2,5,10
a_mc = zeros(1, length(model_orders));
var_mc = zeros(1, length(model_orders));
h_mc = zeros(N/2, length(model_orders));
w_mc = zeros(N/2, length(model_orders));
j = 1;
for i=model_orders
    % Method like above with xcorr
%     a_mc(j) = - Rss(N+i) / Rss(N+i-1);
%     var_mc(j) = Rss(N+i-1) + a1*Rss(N+i);
%     a = [1 a_mc];
%     b = var_mc(j);
%     [h_mc(:,j), w_mc(:,j)] = freqz(b,a,N/2);

    % Aryule method
    [a_mc_yule, var_mc_yule] = aryule(sun_mc, model_orders(j));
    [h_mc_yule(:,j), w_mc_yule(:,j)] = freqz(var_mc_yule^(1/2), a_mc_yule, N/2);
    
    j = j+1;
end


% Plot: og data, normalized, zoomed
f = (1/N):(1/N):((N)/N);
figure('PaperPosition', [0 0 45 15]); 
% og data
subplot(1,3,1); hold on
plot(f, sun_og_psd, 'r', 'DisplayName', 'Periodogram');
for i=1:length(model_orders)
    plot(w_og_yule(:,i)/(2*pi),abs(h_og_yule(:,i)).^2, 'LineWidth', 2, 'DisplayName',strcat('Model PSD AR(', num2str(model_orders(i)),')'));
end
title('Original PSD');
ylabel('Power');
lgd = legend('show'); lgd.Location = 'northeast'; lgd.FontSize = 14;
set(gca,'Xlim',[0,0.5]);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off

% norm data
subplot(1,3,2); hold on
plot(f, sun_mc_psd, 'r', 'DisplayName', 'Periodogram');
for i=1:length(model_orders)
    plot(w_mc_yule(:,i)/(2*pi),abs(h_mc_yule(:,i)).^2, 'LineWidth', 2, 'DisplayName',strcat('Model PSD AR(', num2str(model_orders(i)),')'));
end
title('Normalized PSD');
xlabel('Normalized frequency (Hz)'); 
lgd = legend('show'); lgd.Location = 'east'; lgd.FontSize = 14;
set(gca,'Xlim',[0,0.5]);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off

% zoomed norm
subplot(1,3,3); hold on
plot(f, sun_mc_psd, 'r', 'DisplayName', 'Periodogram');
for i=1:length(model_orders)
    plot(w_mc_yule(:,i)/(2*pi),abs(h_mc_yule(:,i)).^2, 'LineWidth', 2, 'DisplayName',strcat('Model PSD AR(', num2str(model_orders(i)),')'));
end
title('Normalized PSD');
lgd = legend('show'); lgd.Location = 'northwest'; lgd.FontSize = 14;
set(gca,'Xlim',f_zoom);
set(gca, 'Fontsize', 20);
grid on
grid minor
hold off

saveas(gcf, 'CW3_q325', 'epsc');





