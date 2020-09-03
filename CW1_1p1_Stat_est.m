% 1.1 Stat est
clear all; close all;

% Using the MATLAB command rand, generate a 1000-sample vector x = [x[1], x[2],...,x[1000]]T 
% where each sample x[n] is a realisation of a uniform random variable X ? U(0, 1) at time instant n. 
% Plot the result and observe that despite its stochastic nature, x exhibits a 
% degree of uniformity due to its time-invariant statistical properties, since the different
% samples x[n], x[m] have been drawn from the same distribution. 
% Such signals are referred to as statistically stationary.
% The vector x can be considered as a 1000-sample realisation of a stationary stochastic process Xn, whereby Xn ?
% U(0, 1), ?n (all n).


% 1000-sample vector x
x = rand(1000,1);

% Plot the result
% figure('PaperPosition', [0 0 13 10]);
% plot(x); 
% title(['\fontsize{16}Signal of X U(0,1)']); xlabel('n'); ylabel('X[n]'); 
% set(gca, 'Fontsize', 20);



% 1.1.1 Calculate the expected value of X, denoted by m = E{X}, also known as the theoretical mean
% Theo mean of X: integral of x times its probability distribution U(0,1)
m_theo = (1/(1-0)) * ((1^2 - 0^2)/2);        % (1/(b-a)) * (b^2 - a^2)/2 

% Sample mean of X
m_sample = sum(x) / size(x,1);



% 1.1.2 Repeat the analysis for the standard deviation: calculate the theoretical value /sigma = pE{X  E{X}}2 and also its /sighat [5]
% sample estimate from data x using the MATLAB function std which computes the sample standard deviation 

% Theo standard dev (look up derivation this shit long)
    % var = E((X - mew)^2)
    % var = integ((x-mew)^2 * U(0,1) dx)|a,b
    % var = integ((x-mew)^2 * (1/(b-a)) dx)|a,b
    % var = (1/(b-a)) * integ((x-mew)^2 dx)|a,b
    % var = (1/(b-a)) * integ(x^2 - x*mew + mew^2 dx)|a,b
    % var = (1/(b-a)) * [ x^3/3 - mew*x^2/2 + x*mew^2 ]|a,b
    % var = (1/(b-a)) * [ b^3/3 - mew*b^2/2 + b*mew^2  - a^3/3 + mew*a^2/2 - a*mew^2 ]
    % var = (1/(b-a)) * [ (b^3 - a^3)/3 - mew*(b^2 + a^2)/2 + (b - a)*mew^2 ]
    % var = (1/(b-a)) * [ (b^3 - a^3)/3 - mew*(b^2 + a^2)/2 + (b - a)*mew^2 ]
    % etc. 
    % var = E(X^2) - [E(X)]^2
    % var = integ(x^2 * U(0,1) dx) - mew^2
    % var = integ(x^2 * 1/(b-a) dx) - mew^2
    % var = 1/(b-a) * [ x^3/3 ]a,b - mew^2 
    % var = 1/(b-a) * (b^3 - a^3)/3 - ((b+a)/2)^2
    % var = 1/(b-a) * (b^3 - a^3)/3 - (b^2+2ab+a^2)/4
    % etc.
    % var = (b-a)^2 / 12 !!!!!!

std_theo = sqrt( (1-0)^2 / 12 );       % var = E(X^2) - mew^2 = (b-a)^2 / 12

% Sample stand devx
std_sample = std(x);






%% 1.1.3 Generate an ensemble of ten 1000-sample realisations of X, 
% denoted by x1:10, and calculate the sample means mb 1:10 and standard 
% deviations b1:10 for each realisation. Plot these estimates of mean and 
% standard deviation and comment on their bias, by showing how they cluster 
% about their theoretical values.

% Samples x1:10
x_ten = rand(1000,10);

% Sample means and standard dev

std_sample_ten = zeros(1,10); m_sample_ten = zeros(1,10);
for i=(1:10)
    std_sample_ten(i) = std(x_ten(:, i));
    m_sample_ten(i) = mean(x_ten(:, i));
end

% Plot these estimates of mean and standard deviation
% CW1_q113
figure('PaperPosition', [0 0 13 10]);
hold on
plotted_m_theo = ones(size(m_sample_ten,2))*m_theo;
plotted_std_theo = ones(size(std_sample_ten,2))*std_theo;
plot(m_sample_ten, 'r*'); plot(std_sample_ten, 'g*'); 
plot(plotted_m_theo, 'r'); plot(plotted_std_theo, 'g'); 
title(['Statistical estimators for X']); 
xlabel('Realization index n'); ylabel('X[n]'); 
lgd = legend('Mean estimate', 'Stand. dev. estimate'); lgd.Location = 'east';
set(gca, 'Fontsize', 20);
txt_m = ['Var of mean estimator: ', num2str((std(m_sample_ten))^2)];
txt_std = ['Var of std estimator: ', num2str((std(std_sample_ten))^2)];
text(2,0.35,txt_m,'FontSize',15); text(2,0.32,txt_std,'FontSize',15);
hold off
saveas(gcf, 'CW1_q113', 'epsc');





%% 1.1.4 Approximate the pdf of X by showing in the same plot the histogram of x (see hist),
% normalised by the number of samples considered, and the theoretical pdf. 
% Comment upon the result, in particular, on whether the estimate appears 
% to converge as the number of generated samples increases, and how the number of
% histogram bins considered affects the analysis.

% Histogram of x: CW1_q114
figure('PaperPosition', [0 0 40 10]);
% Varying sample size
sgtitle(['Histogram and fit of X U(0,1)']); 
binsize = 100;
plotcount = 1;
for sample_num = [100,1000,2000]
    subplot(1, 3,plotcount); 
    x = rand(sample_num,1);
    unif_pdf = (ones(size(x,1)))./binsize;
    hold on
    histogram(x, binsize, 'Normalization','probability');
    plot(x, unif_pdf, 'r');
    hold off
    title(['Sample size=', num2str(sample_num)]); 
    xlabel('binned X[n]'); ylabel('Relative counts'); 
    lgd = legend('Histogram X', 'Theoretical pdf'); lgd.FontSize = 15;
    set(gca, 'Fontsize', 20);
    plotcount = plotcount+1;
end
saveas(gcf, 'CW1_q114_unif', 'epsc');

%% Varying histograms
figure('PaperPosition', [0 0 40 10]);
sgtitle('Histogram and fit of X U(0,1)'); 
x = rand(1000,1);
plotcount = 1;
for binsize = [20,100,500]
    unif_pdf = (ones(size(x,1)))./binsize;
    subplot(1, 3,plotcount); 
    hold on
    histogram(x, binsize, 'Normalization','probability');
    plot(x, unif_pdf, 'r');
    hold off
    title(['Binsize=', num2str(binsize)]); 
    xlabel('binned X[n]'); ylabel('Relative counts'); 
    lgd = legend('Histogram X', 'Theoretical pdf'); lgd.FontSize = 15; 
    if plotcount ==1
        lgd.Location = 'west';
    end
    set(gca, 'Fontsize', 20);
    plotcount = plotcount+1;
end
saveas(gcf, 'CW1_q114_bins', 'epsc');



%% 1.1.5 REPEAT 1-4 LOL HAHAH FUCK YOU use randn and get N(0,1)

% 1000-sample vector x
y = randn(1000,1);

% Plot the result
% figure;
% plot(y); 
% title(['\fontsize{16}Signal of X N(0,1)']); xlabel('n'); ylabel('X[n]'); 


% 1.1.5.1 
% Theo mean of X N(0,1)
m_ytheo =  0;        % (1/(b-a)) * (b^2 - a^2)/2 

% Sample mean of X
m_ysample = mean(y);



% 1.1.5.2 std
std_ytheo = 1;       % var = E(X^2) - mew^2 = (b-a)^2 / 12;

% Sample stand devx
std_ysample = std(y);


% 1.1.5.3 generate  x1:10
% Samples x1:10
y_ten = randn(1000,10);


%% Sample means and standard dev

std_ysample_ten = zeros(1,10); m_ysample_ten = zeros(1,10);

for i=(1:10)
    std_ysample_ten(i) = std(y_ten(:, i));
    m_ysample_ten(i) = mean(y_ten(:, i));
end

% Plot these estimates of mean and standard deviation
figure('PaperPosition', [0 0 12 10]);
hold on
plotted_m_ytheo = ones(size(m_ysample_ten,2))*m_ytheo;
plotted_std_ytheo = ones(size(std_ysample_ten,2))*std_ytheo;
plot(m_ysample_ten, 'r*'); plot(std_ysample_ten, 'g*'); 
plot(plotted_m_ytheo, 'r'); plot(plotted_std_ytheo, 'g'); 
title(['Statistical estimators for X N(0,1)']); 
xlabel('Realization index n'); ylabel('X[n]'); 
lgd = legend('Mean estimate', 'Stand. dev. estimate'); lgd.Location = 'east';
txt_m = ['Var of mean est: ', num2str((std(m_sample_ten))^2)];
txt_std = ['Var of std est: ', num2str((std(std_sample_ten))^2)];
text(2,0.3,txt_m,'FontSize',15); text(2,0.15,txt_std,'FontSize',15);
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW1_q113b', 'epsc');



%% 1.1.5.4 Hist
% Histogram of x N(0,1)
figure('PaperPosition', [0 0 35 10]);
binsize_y = 100;
% histogram(y, binsize_y, 'Normalization','probability'); 
% norm_pd = fitdist((y./binsize_y),'normal');
% dummy_y = 0.01:0.01:1; norm_curve = pdf(norm_pd, dummy_y);
% plot(norm_curve);
plot_count = 1;
for sample_size = [100,1000,2000]
    
    hold on
    subplot(1,3, plot_count);
    y = randn(sample_size, 1);
    histfit(y, binsize_y, 'normal');
    title(['Sample num=', num2str(sample_size)]); 
    xlabel('binned X[n]'); ylabel('Relative counts'); 
    lgd = legend('Histogram X', 'Theoretical pdf'); lgd.FontSize = 15;
    yt = get(gca, 'YTick');
    set(gca, 'YTick', yt, 'YTickLabel', yt/numel(y));
    set(gca, 'Fontsize', 20);
    hold off
    plot_count=plot_count+1;

end
saveas(gcf, 'CW1_q114_binsb', 'epsc');


