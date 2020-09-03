% ASP CW 3.5
clear all 
close all

%% 3.5.1

% Turn ECG to RRI
load('xRRI_trial1.mat');
load('xRRI_trial2.mat');
load('xRRI_trial3.mat');

% Normalize (zscore)
xRRI_trial1 = detrend(xRRI_trial1);
xRRI_trial2 = detrend(xRRI_trial2);
xRRI_trial3 = detrend(xRRI_trial3);

N = size(xRRI_trial1, 2);

st_psd = zeros(N, 3);
st_psd(:,1) = pgm(xRRI_trial1);
st_psd(:,2) = pgm(xRRI_trial2);
st_psd(:,3) = pgm(xRRI_trial3);

f = 0:(1/N):((N/2-1)/N);

figure('PaperPosition', [0 0 14 10]); 
plot(f,st_psd(1:N/2,1)); 
title(['RRI PSD Trial 1']);
xlabel('Normalized Frequency (Hz)'); ylabel('Power');
set(gca, 'Fontsize', 20); 
saveas(gcf, 'CW3_q351_1', 'epsc');

figure('PaperPosition', [0 0 14 10]); 
plot(f,st_psd(1:N/2,2)); 
title(['RRI PSD Trial 2']);
xlabel('Normalized Frequency (Hz)'); ylabel('Power');
set(gca, 'Fontsize', 20); 
saveas(gcf, 'CW3_q351_2', 'epsc');

figure('PaperPosition', [0 0 14 10]); 
plot(f,st_psd(1:N/2,3));
title(['RRI PSD Trial 3']);
xlabel('Normalized Frequency (Hz)'); ylabel('Power');
set(gca, 'Fontsize', 20); 
saveas(gcf, 'CW3_q351_3', 'epsc');


% Ave PSD's
N = size(xRRI_trial1, 2);
w_lengths = [N/20, N/10, N/5];     % window lengths (s)
% psd_win = zeros(N, 3);
% psd_ave = zeros(size(psd_win));
plot_num = 1;
figure('PaperPosition', [0 0 28 10]);
for curr_win = w_lengths
%     psd_mean = zeros(N/curr_win,3);
    psd_wind = zeros(curr_win, N/curr_win,3);
    for div = 1:N/curr_win
        i1 = (div-1)*curr_win+1;
        i2 = div*curr_win;
%         psd_win(i1:i2, 1) = pgm(xECG_trial1(i1:i2));
%         psd_win(i1:i2, 2) = pgm(xECG_trial2(i1:i2));
%         psd_win(i1:i2, 3) = pgm(xECG_trial3(i1:i2));
        
        psd_wind(1:curr_win, div,1) = pgm(xRRI_trial1(i1:i2));
        psd_wind(1:curr_win, div,2) = pgm(xRRI_trial2(i1:i2));
        psd_wind(1:curr_win, div,3) = pgm(xRRI_trial3(i1:i2));
        
%         psd_mean(i,:) = mean(psd_win(i1:i2, :));
    end
    psd_mean = mean(psd_wind,2);
    
    f = 0:(1/curr_win):((curr_win/2-1)/curr_win);

    subplot(1,3,plot_num); hold on
    plot(f, (psd_mean(1:floor(curr_win/2),:,1))); 
    plot(f, (psd_mean(1:floor(curr_win/2),:,2))); 
    plot(f, (psd_mean(1:floor(curr_win/2),:,3)));
    title(['Ave PSD win=', num2str(curr_win)]);
    xlabel('Norm Freq (Hz)'); ylabel('Power');
    legend('Trial 1', 'Trial 2', 'Trial 3');
    set(gca, 'Fontsize', 20); 
    hold off    
    
    plot_num = plot_num+1;
end
saveas(gcf, 'CW3_q351_ave', 'epsc');

