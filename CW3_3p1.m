% CW 3.0 and 3.1
clear all; close all; 

% 3.0.1: Periodogram of WGN

N = [128, 256, 512];
figure('PaperPosition', [0 0 30 8]); hold on
plot_count = 1;
for n = N
    wgn = randn(1, n);
    periodogram = pgm(wgn);
    f = 0:(1/n):((n-1)/n);
    
    subplot(1,length(N), plot_count); 
    plot(f, 10*log10(periodogram));
    title(['PSD of WGN, N=', num2str(n)]);
    xlabel('Normalized frequency (Hz)'); ylabel('Power (dB)');
    set(gca, 'Fontsize', 20);
    
    plot_count = plot_count +1;
end
hold off
saveas(gcf, 'CW3_q301', 'epsc');


%% 3.1.1 Filter PSD

N = 256;
b = 0.2*[1 1 1 1 1];
a = 1;
f = 0:(1/N):((N-1)/N);

wgn = randn(1,N);
pgm_og = pgm(wgn);
pgm_filt = filter(b, a, pgm_og);

figure; hold on
plot(f, 10*log10(pgm_og),'LineWidth',2); plot(f, 10*log10(pgm_filt),'LineWidth',2);
title(['PSD of WGN, N=', num2str(N)]);
xlabel('Normalized frequency (Hz)'); ylabel('Power (dB)');
legend('Original PSD', 'Smooth PSD');
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW3_q311', 'epsc');

var_og_pgm = var(pgm_og);      % variance of pgm 
var_smooth_pgm = var(pgm_filt);
    

%% 3.1.2 Subdivide & PSD
% [sample id, data]
N = 1024;
N2 = 128;
wgn = randn(1,N);
wgn_div = reshape(wgn, [8, N2]);       % 8 samples, 128 data length

% get psd
pgm_div = zeros(size(wgn_div));         
for div = 1:size(wgn_div,1)
    pgm_div(div,:) = pgm(wgn_div(div,:));
end

var_pgm = var(pgm_div');        % transpose so that the var of each segment found

% plot variance
figure; 
bar(var_pgm, 'LineWidth', 2); 
xlabel('Sample segment'); ylabel('Periodogram variance (dB)');
title('Variance of WGN Periodogram');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW3_q312', 'epsc'); 


%% 3.1.3
pgm_ave_est = mean(pgm_div,1);

% display average of all segments
f = 0:(1/N2):((N2-1)/N2);

figure; hold on 
for i=1:8
    plot(f, 10*log10(pgm_div(i,:)), 'b');
end
plot(f, 10*log10(pgm_ave_est), 'r', 'LineWidth', 2);
title('Averaged Periodogram'); 
xlabel('Normalized frequency (Hz)'); ylabel('Power (dB)');
legend('segment 1', 'segment 2','segment 3', 'segment 4','segment 5','segment 6', 'segment 7', 'segment 8', 'segment ave');
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW3_q313', 'epsc');





