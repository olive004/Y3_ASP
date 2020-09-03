% ASP CW 4.4 
clear all;
close all;

% rng(0);     % set random num gen seed


%% 4.4.1 Block diagram

% Synthesis
% wgn
N = 4000;
x = randn(N, 1);

% AR model
a = [1, 0.9, 0.2];
b = [1];
x = filter(b, a, x);
% x = zscore(x);
% x = x./sqrt(sum(a.*a));
x_decorr = [0; x(1:N-1)];

% Analysis
N_w = length(a);
filt_order = N_w-1;

mu = 0.01;
% [x_hat, e_evol, w_evol] = lms(x_decorr, x, mu, filt_order);
[x_hat, e_evol, w_evol] = lms_adap(x, mu, filt_order);

% plot error
sq_e = e_evol.^2;
figure('PaperPosition', [0 0 25 7]);
subplot(1,2,1); 
plot(sq_e, 'Linewidth', 2);
title(['LMS Error^2']);
xlabel('Iteration n'); ylabel('error');
set(gca, 'Fontsize', 20);

% plot weight evol
subplot(1,2,2); 
plot(w_evol');
title(['Weight evolution, mu=', num2str(mu)]);
xlabel('Iteration n'); ylabel('weights');
legend('a1', 'a2');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q441', 'epsc');

disp(w_evol(:,end));


%% 4.4.2 No gear shifting

lrs = [0.001, 0.01, 0.1, 0.2];

% plot weight evol
figure('PaperPosition', [0 0 25 22]);
plot_num = 1;
for mu = lrs 
    [x_hat, e_evol, w_evol] = lms_adap(x, mu, filt_order);
    a1(plot_num) = w_evol(1,end);
    a2(plot_num) = w_evol(2,end);
    
    subplot(2,2,plot_num); 
    plot(w_evol');
    title(['Weight evolution, mu=', num2str(mu)]);
    xlabel('Iteration n'); ylabel('weights');
    legend('a1', 'a2');
    set(gca, 'Fontsize', 20);
    
    plot_num = plot_num+1;
    
%     subplot(2,4,plot_num); 
%     hold on
%     plot(x); plot(x_hat);
%     title(['Adapted signal']);
%     xlabel('Iteration n'); ylabel('weights');
%     legend('x', 'x_{hat}');
%     set(gca, 'Fontsize', 20);
% 
%     plot_num = plot_num+1;
end
disp(a1); disp(a2);
saveas(gcf, 'CW4_q442', 'epsc');

