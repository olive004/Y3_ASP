% Q 1.3
clear all; close all;

% 1.3.1
N = 1000;
n_bins = 100; 
smoothing_factor = 5;
v = randn(1,N);

figure;
pdf(v, n_bins, smoothing_factor);
title(['Estimated pdf v at N=', num2str(N)]); 
xlabel('v[n]'); ylabel('Relative counts');
lgd = legend("Histogram", "Estimated pdf"); lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW1_q131_hist', 'epsc');


% Code from 1.2 needed for this part:
M = 100;
NN = [100 1000 10000];


figure('PaperPosition', [0 0 30 8]); hold on
plot_count = 1;
smoothing_factor = 7;
sgtitle(['Estimated pdf x3']);
for N = NN
    x3 = rp3(M, N);
    hAx = subplot(1,3,plot_count);
    pdf(x3, n_bins, smoothing_factor);
    
    xlabel('X_3[n]'); ylabel('Relative counts'); 
    ylim(hAx, [0, 0.5]);
    title(['N=', num2str(N)]); 
    set(gca, 'Fontsize', 20);
    plot_count = plot_count+1;
    lgd = legend("Histogram", "Estimated pdf"); lgd.FontSize = 15; lgd.Location = 'southwest';
end

hold off
saveas(gcf, 'CW1_q132_HistsN', 'epsc');

% 1.3.3 The big boi
% PDF of one of the nonstationary function rp1
N = 1000;
M =1;
n_bins = 20;
x1 = rp1(M, N);

figure; 
smoothing_factor = 3;
pdf(x1, n_bins, smoothing_factor);
title(['Estimated pdf of nonstationary rp1, N=', num2str(N)]); 
xlabel('X_1[n]'); ylabel('Relative counts'); 
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW1_q133_nonstat', 'epsc');

figure; plot(x1);
title(['Signal rp2, N=', num2str(N)]); 
xlabel('Index n'); ylabel('X_1[n]'); 
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW1_q133_signal', 'epsc');


x_piecewise = [rand(M, N), (rand(M, N) +1)];
figure; 
pdf(x_piecewise, n_bins, smoothing_factor);
title(['Nonstationary signal with mean=[0,1], N=', num2str(N)]); 
xlabel('X_1[n]'); ylabel('Relative counts'); 
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW1_q133_piecewise', 'epsc');


% for N = NN
% %     x1 = rp1(M, N);
% %     x2 = rp2(M, N);
%     x3 = rp3(M, N);
% 
%     %%%%%%%%%%%%%%%%%%%
% %     ONLY RP3 IS ERG AND STAT!!!!!
%     % Q 1.3.2   For those processes in Part 1.2 (a,b,c) that are stationary and
%     % ergodic, run your MATLAB code to approximate the pdf for N 2 {100, 1000, 
%     % 10000}. Compare the estimated densities with the theoretical pdf using the 
%     % MATLAB subplot function, and comment on the results as data length N increases.
% 
%     % Only x2 and x3 are ergodic and stationary
% 
% %     [pdf2_cent, pdf2_counts] = pdf_c(x2, n_bins);
%     [pdf3_cent, pdf3_counts] = pdf_c(x3, n_bins);
% 
%     % Compare est & theo pdfs
%     % x2 has pdf = U[0,1]^2 + 0.5 U[0,1] = (integ(x)|0,1)^2 + 0.5*integ(x)|0,1   
%     % pdf2 = (x^2/2)(x^2/2) + 0.5((x^2/2)
%     % pdf2 = (x^4/4) + (x^2/4) = (x^2/4)(x^2 + 1)  |0,1
% %     pdf2_theo = fitdist(pdf2_cent,'Normal', 'lower',-1,'upper',2);
% %     pdf2_theo_sampled = pdf(pdf2_theo, pdf2_cent);
% 
%     % x3 has pdf = 3U[-0.5,0.5] + 0.5 
%     pdf3_theo = makedist('Uniform', 'lower',-1,'upper',2);
%     pdf3_theo_sampled = pdf(pdf3_theo, pdf3_cent);
% 
%     % Plotting estimated and theo pdfs
%     % pdf 2
% %     figure; subplot(2,1,1); histogram(x2, 'Normalization', 'pdf'); 
% %     title(['Estimated pdf x2 at N=', num2str(N)]); xlabel('x2 values'); ylabel('counts');
% %     subplot(2,1,2); plot(pdf2_cent, pdf2_theo_sampled);
% %     title('Theoretical pdf x2'); xlabel('x2 values'); ylabel('counts');
% 
%     % pdf 3
%     figure; 
%     hAx(1) = subplot(2,1,1); histogram(x3, 'Normalization', 'pdf'); 
%     title(['Estimated pdf x3 at N=', num2str(N)]); xlabel('x3 values'); ylabel('counts');
%     hAx(2) = subplot(2,1,2); plot(pdf3_cent, pdf3_theo_sampled);
%     title('Theoretical pdf x3'); xlabel('x3 values'); ylabel('counts');
%     ylim(hAx, [0, 0.5])
% 
% 
%     
% end






