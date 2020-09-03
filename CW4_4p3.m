% ASP CW4_4p3 Gear shifting
clear all;
close all;

rng(0);     % set random num gen seed


% 4.1 & 4.2
b = [1, 2, 3, 2, 1];
a = [1];
N_w = length(b);         % num coefficients
filt_order = N_w + 1;
N = 1000; %filt_order;
x = randn(N,1);

% unknown system
y = filter(b, a, x);
y = zscore(y);

std = 0.1;
nu = std*randn(N,1);       % wgn 2

z = y + nu;


initial_lrs = [0.01, 0.02, 0.2, 1];
i = 1;
for mu = initial_lrs
    [y_hat, e, weight_evol] = lms(x, z, mu, filt_order);

    % plot error
    sq_e = e.^2;
    figure('PaperPosition', [0 0 25 7]);
    subplot(1,2,1); 
    plot(sq_e, 'Linewidth', 2);
    title(['LMS Error^2, mu=', num2str(mu)]);
    xlabel('Iteration n'); ylabel('error');
    set(gca, 'Fontsize', 20);
    
    % plot weight evol
    subplot(1,2,2); 
    plot(weight_evol');
    title(['Weight evolution, mu=', num2str(mu)]);
    xlabel('Iteration n'); ylabel('weights');
    set(gca, 'Fontsize', 20);
    plot_name = ['CW4_q431_', num2str(i)];
    saveas(gcf, plot_name, 'epsc');
    
    disp(i); disp(weight_evol(:,end));
    
    i = i+1;
end

