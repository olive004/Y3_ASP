% ASP CW 4.1 
clear all
close all


N = 1000;
x = randn(N,1);

% unknown system
b = [1, 2, 3, 2, 1];
a = [1];
y = filter(b, a, x);

y = zscore(y);

std = 0.1;
nu = std*randn(N,1);       % wgn 2
z = y + nu;

snr1 = snr(z, nu);


%% 4.1.1

Nw = length(b);             % number of w_opt coeffs
i1 = N; i2 = Nw-1;
Rxx = xcorr(x,'unbiased');
Rxx = toeplitz(Rxx(i1:i1+i2));
pzx = xcorr(z,x, 'unbiased');
figure; plot(pzx);
pzx = pzx(i1:i1+i2);

w_opt = Rxx\pzx; 

% if (sum(find(w_opt)) <= 2) || (w_opt(1)<0)      % w_opt(1)<0 || w_opt(2)<0 || w_opt(3)<0 || w_opt(4)<0
%     y_est = filter(-w_opt, a, x);
% else 
%     y_est = filter(w_opt, a, x);
% end
y_est = filter(w_opt, a, x);
y_est = zscore(y_est);
y = zscore(y);

figure('PaperPosition', [0 0 12 8]); hold on
plot(y);
plot(y_est);
title('Unknown and estimated y');
xlabel('n'); ylabel('strength');
lgd = legend('y[n]', 'estimated y[n]'); lgd.Location = 'south';
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q411', 'epsc');


%% 4.1.2 

% Get PSDs
window_size = N/10;
psd_y = get_ave_psd(y,window_size);

figure('PaperPosition', [0 0 35 25]);
stds = linspace(0.1, 10, 6);
% and error
mse = zeros(1, length(stds));
snr_z = zeros(1,length(stds));
i1 = N; i2 = Nw-1;
plot_num = 1;
for std = stds
    nu = std*randn(N,1);       % wgn 2
    z = y + nu;

    snr_z(plot_num) = snr(z, nu);
    
    Rxx = xcorr(x,'unbiased');
    Rxx = toeplitz(Rxx(i1:i1+i2));
    pzx = xcorr(z,x, 'unbiased');
    pzx = pzx(i1:i1+i2);

    w_opt = Rxx\pzx; 

    y_est = filter(w_opt, a, x);
    y_est = zscore(y_est);
    psd_est = get_ave_psd(y_est,window_size);
    
    mse(plot_num) = sum((y_est - y).^2);
    
    
    f = 0:(1/window_size):(window_size/2-1)/window_size;
    subplot(length(stds)/3,length(stds)/2, plot_num); hold on
    plot(f, psd_y(1:window_size/2));
    plot(f, psd_est(1:window_size/2));
    title(['PSD y, std=', num2str(std)]);
    xlabel('Norm freq (Hz)'); ylabel('Power (dB)');
    legend('y[n]', 'estimated y[n]');
    set(gca, 'Fontsize', 20);

    plot_num = plot_num +1;
end
saveas(gcf, 'CW4_q412', 'epsc');

% normalize MSE
mse = mse./max(mse);

figure;
stem(stds, mse);
title('MSE vs std'); 
xlabel('Stand. dev.'); ylabel('Normalized MSE');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q412_mse', 'epsc');
