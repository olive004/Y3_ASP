% CW4 4.2 LMS
clear all
close all

%% 4.2.1

% FOLLOWING PARTLY FROM 4.1
b = [1, 2, 3, 2, 1];
a = [1];
Nw = length(b);         % num coefficients
filt_order = Nw + 1;
N = 1000; %filt_order;
x = randn(N,1);

% unknown system
y = filter(b, a, x);
y = zscore(y);

std = 0.1;
nu = std*randn(N,1);       % wgn 2

z = y + nu;



% apply LMS function
mu = 1;             % adaptive gain
[y_hat, e, weight_evol] = lms(x, z, mu, filt_order);
y_hat = zscore(y_hat);

% Comparison y to yhat
sq_e = e.^2;

figure; hold on
plot(y); 
plot(y_hat);
title(['LMS estimated signal, filter order=', num2str(filt_order)]);
xlabel('n'); ylabel('Signal strength');
legend('y', 'y_{hat}');
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW4_q421', 'epsc');

% plot error
figure; plot(sq_e, 'Linewidth', 2);
title(['LMS Error^2, mu=', num2str(mu)]);
xlabel('Iteration n'); ylabel('error');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q421_e', 'epsc');

% plot weights
figure; plot((N-50:N), weight_evol(:,N-50:N)');
title(['Weight evolution, mu=', num2str(mu)]);
xlabel('Iteration n'); ylabel('weights');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q421_w', 'epsc');


%% 4.2.2 Lower lr

lrs = [0.005, 0.01, 0.1, 0.25];
i = 1;
for mu = lrs
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
    plot_name = ['CW4_q422_', num2str(i)];
    saveas(gcf, plot_name, 'epsc');
    
    disp(i); disp(weight_evol(:,end));
    
    i = i+1;
end


