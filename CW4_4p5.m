% ASP CW 4: 4.5
clear all;
close all;


%% 4.5

% come in 2 identical columns, drop 1
e = audioread('e.wav');
a = audioread('a.wav');
s = audioread('s.wav');
t = audioread('t.wav');
x = audioread('x.wav');

info = [audioinfo('e.wav'); %get the wav information
        audioinfo('a.wav')
        audioinfo('s.wav')
        audioinfo('t.wav')
        audioinfo('x.wav')];



% figure; hold on
% plot(x); plot(movmean(abs(x), 400));
% b = (1/window_size)*ones(1,window_size);
% a = 1;
% t_f = filter(b,a,x);
% plot(t_f);

% Delete silence
window_size = 300;
old_length = [length(e), length(a), length(s), length(t), length(x)];

e = del_silence(e, window_size, 1);
a = del_silence(a, window_size, 1);
s = del_silence(s, window_size, 1);
t = del_silence(t, window_size, 4500);
x = del_silence(x, window_size, 1);

new_length = [length(e), length(a), length(s), length(t), length(x)];
new_duration = [(info(1).Duration * new_length(1)) / old_length(1),
    (info(2).Duration * new_length(2)) / old_length(2),
    (info(3).Duration * new_length(3)) / old_length(3),
    (info(4).Duration * new_length(4)) / old_length(4),
    (info(5).Duration * new_length(5)) / old_length(5)];

%% 4.5.1

% Adaptation
fs = 44100;
N = 1000;

% Resample audio
e_samp = resample(e, fs, info(1).SampleRate); 
a_samp = resample(a, fs, info(2).SampleRate); 
s_samp = resample(s, fs, info(3).SampleRate);
t_samp = resample(t, fs, info(4).SampleRate); 
x_samp = resample(x, fs, info(5).SampleRate); 

% Use certain sample
to_end=1;
if to_end
    e_samp = e_samp(1:end);
    a_samp = a_samp(1:end);
    s_samp = s_samp(1:end);
    t_samp = t_samp(1:end);
    x_samp = x_samp(1:end);
else
    e_samp = e_samp(1:N);
    a_samp = a_samp(1:N);
    s_samp = s_samp(1:N);
    t_samp = t_samp(1:N);
    x_samp = x_samp(1:N);
end

e_samp= zscore(e_samp); a_samp= zscore(a_samp); s_samp= zscore(s_samp);
t_samp= zscore(t_samp); x_samp= zscore(x_samp);



%% 'e' 1

% Test model orders
Nw_max = 10;
filt_orders = 1:Nw_max;
lr = 0.01;

z_all_w = zeros(Nw_max, length(filt_orders));
time_e = new_duration(1)*(0:1/length(e_samp):(length(e_samp)-1)/length(e_samp));
figure('PaperPosition', [0 0 24 10]);
c={'r','g','b','k','c','m'};
w = 1;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(e_samp, lr, filt_order, 0);
    z_all_w(w, :) = [z_w(:, end);  zeros(1, (length(filt_orders)-filt_order))']';
    
    % Split into 2 subplots
    if w>length(filt_orders)/2 
        subplot(1,2,2); 
        hold on
        p = plot(z_w','color',c{w - (length(c)-1)*floor(w/length(c))});
        for j = 2:w
            set(get(get(p(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        ylim([0 1]);
        title(['Weights "e", N_w=', num2str(length(filt_orders)/2+1),':',num2str(filt_order)]);
        xlabel('Iterations N'); ylabel('Weight');
        set(gca, 'Fontsize', 20);
        lgd = legend('w(6)', 'w(7)', 'w(8)', 'w(9)', 'w(10)'); lgd.Location = 'northwest'; lgd.FontSize = 15;
        hold off
    else
        subplot(1,2,1); 
        hold on
        p = plot(z_w','color',c{w - (length(c)-1)*floor(w/length(c))}); %'LineWidth', 2,
        for j = 2:w
            set(get(get(p(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        ylim([0 1]);
        title(['Weights "e", N_w=', num2str(1),':',num2str(filt_order)]);
        xlabel('Iterations N'); ylabel('Weight');
        set(gca, 'Fontsize', 20');
        lgd = legend('w(1)', 'w(2)', 'w(3)', 'w(4)', 'w(5)'); lgd.Location = 'northwest'; lgd.FontSize = 15;
        hold off
    end
    w= w+1;
end
disp(z_all_w);
saveas(gcf, 'CW4_q451_nw', 'epsc');


%% 'e' 1

% Vary Learning rates
lrs = [0.01, 0.05, 0.1, 0.2];
filt_order = 2;
z_all_w = zeros(length(lrs), filt_order);
w = 1;
c={'m','g','c','k'};
figure('PaperPosition', [0 0 13 10]); hold on
for lr = lrs
    
    [z_hat, z_e, z_w] = lms_adap(e_samp, lr, filt_order, 0);
    
%     plot(e_hat); plot(e);
%     title(['Audio "e", lr=', num2str(lr)]); 
%     xlabel('Time'); ylabel('Amplitude');
%     legend('estimated e', 'e');

%     plot(e_e); 
%     title('Error'); xlabel('Iterations N'); ylabel('Error');
    p = plot(z_w', 'color',c{w});
    set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    z_all_w(w, :) = z_w(:, end);
    
    w = w+1;
end
title('Weight evolution "e"'); xlabel('Iterations N'); ylabel('Weight');
lgd = legend(num2str(lrs(1)), num2str(lrs(2)),num2str(lrs(3)),num2str(lrs(4)));
lgd.Location = 'northwest'; lgd.FontSize = 15;
title(lgd, 'Learning rate');
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW4_q451', 'epsc');

disp(z_all_w);
lrs = [lrs', lrs'];
z_all_w = (1-lrs).*z_all_w;
disp(mean(z_all_w,1));


%% 'e' 1

% With gear shift turned on
lr = 0.01;
filt_order = 2;
[z_hat, z_e, z_w] = lms_adap(e_samp, lr, filt_order, 1);
figure('PaperPosition', [0 0 13 10]);
plot(z_w', 'LineWidth', 2);
title('Gear shift weight, lr=0.006'); xlabel('Iterations N'); ylabel('Weight');
lgd = legend('w(1)', 'w(2)'); lgd.Location = 'northwest'; lgd.FontSize = 15;
title(lgd, 'Weight');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q451_w', 'epsc');


%% 'e' 1

% Show fit of signal with optimal coeffs

figure('PaperPosition', [0 0 20 10]); subplot(1,2,1);
hold on
plot(time_e, e_samp); plot(time_e,z_hat); 
title(['Adapted signal "e", N_w=', num2str(filt_order)]);
xlabel('time (s)'); ylabel('Signal amplitude');
set(gca, 'Fontsize', 20);
lgd = legend('e', 'e_{hat}'); lgd.FontSize = 15;
hold off

subplot(1,2,2); 
plot(z_e);
title(['Error "e"']);
xlabel('Iterations N'); ylabel('Error');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW4_q451_e', 'epsc');









%% 's' etc

% Test model orders 
Nw_max = 10;
filt_orders = 1:Nw_max;
lr = 0.001;

letter = x_samp;

figure('PaperPosition', [0 0 24 10]);
c={'r','g','b','k','c','m'};
z_all_w = zeros(Nw_max, length(filt_orders));
w = 1;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    z_all_w(w, :) = [z_w(:, end);  zeros(1, (length(filt_orders)-filt_order))']';
    
    % Split into 2 subplots
    if w>length(filt_orders)/2 
        subplot(1,2,2); 
        hold on
        p = plot(z_w','color',c{w - (length(c)-1)*floor(w/length(c))});
        for j = 2:w
            set(get(get(p(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        ylim([0 1]);
        title(['Weights, N_w=', num2str(length(filt_orders)/2+1),':',num2str(filt_order)]);
        xlabel('Iterations N'); ylabel('Weight');
        set(gca, 'Fontsize', 20);
        lgd = legend('w(6)', 'w(7)', 'w(8)', 'w(9)', 'w(10)'); lgd.Location = 'northwest'; lgd.FontSize = 15;
        hold off
    else
        subplot(1,2,1); 
        hold on
        p = plot(z_w','color',c{w - (length(c)-1)*floor(w/length(c))}); %'LineWidth', 2,
        for j = 2:w
            set(get(get(p(j),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        ylim([0 1]);
        title(['Weights, N_w=', num2str(1),':',num2str(filt_order)]);
        xlabel('Iterations N'); ylabel('Weight');
        set(gca, 'Fontsize', 20');
        lgd = legend('w(1)', 'w(2)', 'w(3)', 'w(4)', 'w(5)'); lgd.Location = 'northwest'; lgd.FontSize = 15;
        hold off
    end
    
    w= w+1;
end
disp('weights by filter order');
disp(z_all_w);

% Vary Learning rates
lrs = [0.0001, 0.0005, 0.001, 0.005];
filt_order = 2;
z_all_w = zeros(length(lrs), filt_order);
w = 1;
c={'m','g','c','k','r','b'};

figure('PaperPosition', [0 0 13 10]); hold on
for lr = lrs
    [z_hat, z_e, z_w] = lms_adap(letter, lr, filt_order, 0);
    p = plot(z_w', 'color',c{w - (length(c)-1)*floor(w/length(c))});
    set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    z_all_w(w, :) = z_w(:, end);
    
    w = w+1;
end
title('Weight evolution'); xlabel('Iterations N'); ylabel('Weight');
lgd = legend(num2str(lrs(1)), num2str(lrs(2)),num2str(lrs(3)),num2str(lrs(4)));
lgd.Location = 'northwest'; lgd.FontSize = 15;
title(lgd, 'Learning rate');
set(gca, 'Fontsize', 20);
hold off


% With gear shift turned on
lr = 0.0005;
filt_order = 2;
[z_hat, z_e, z_w] = lms_adap(letter, lr, filt_order, 1);
figure('PaperPosition', [0 0 13 10]);
plot(z_w', 'LineWidth', 2);
title('Gear shift weight'); xlabel('Iterations N'); ylabel('Weight');
lgd = legend('w(1)', 'w(2)'); lgd.Location = 'northwest'; lgd.FontSize = 15;
title(lgd, 'Weight');
set(gca, 'Fontsize', 20);
saveas(gcf,'CW4_q451_a_gs','epsc');
disp('gs weights')
disp(z_w(:,end));







%% 4.5.2 Prove AR processes

figure('PaperPosition', [0 0 10 18]);
subplot(5,1,1); 
plot(xcorr(e_samp, 'unbiased')); 
title('ACF "e"'); xlabel('Time lag \tau'); ylabel('Correlation'); 
set(gca, 'Fontsize', 20);
subplot(5,1,2); 
plot(xcorr(a_samp, 'unbiased')); 
title('ACF "a"'); xlabel('Time lag \tau'); 
set(gca, 'Fontsize', 20);
subplot(5,1,3); 
plot(xcorr(s_samp, 'unbiased')); 
title('ACF "s"'); xlabel('Time lag \tau'); ylabel('Correlation'); 
set(gca, 'Fontsize', 20);
subplot(5,1,4); 
plot(xcorr(t_samp, 'unbiased')); 
title('ACF "t"'); xlabel('Time lag \tau'); 
set(gca, 'Fontsize', 20);
subplot(5,1,5); 
plot(xcorr(x_samp, 'unbiased')); 
title('ACF "x"'); xlabel('Time lag \tau'); ylabel('Correlation'); 
set(gca, 'Fontsize', 20);
saveas(gcf,'CW4_q452_acf','epsc');




%% 4.5.2 R_p and Optimal filter length

N = 1000;
Nw_max = 10;
filt_orders = 1:Nw_max;
lr = 0.001;

R_p = zeros(5, length(filt_orders));
errors = zeros(5, length(filt_orders));

w = 1;
letter = e_samp; lr = 0.01;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(1,w) = var(z_e);
    
    var_x = var(letter);
    var_e = var(z_e);
    R_p(1, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = a_samp; lr = 0.001;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(2,w) = var(z_e);
    
    var_x = var(letter);
    var_e = var(z_e);
    R_p(2, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = s_samp; lr = 0.005;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(3,w) = var(z_e);

    var_x = var(letter);
    var_e = var(z_e);
    R_p(3, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = t_samp; lr = 0.001;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(4,w) = var(z_e);

    var_x = var(letter);
    var_e = var(z_e);
    R_p(4, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = x_samp; lr = 0.001;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(5,w) = var(z_e);

    var_x = var(letter);
    var_e = var(z_e);
    R_p(5, w) = 10*log10(var_x / var_e);
    
    w= w+1;
end

% Normalize
for w=1:5
    R_p(w, :) = R_p(w, :)./max(R_p(w, :));
%     errors(i, :) = errors(i, :)./max(errors(i, :));
end

figure('PaperPosition', [0 0 30 10]);
subplot(1,2,1);
plot(R_p'); 
title('R_p by filter order'); xlabel('Filter order'); ylabel('R_p'); xlim([-2 length(filt_orders)]);
lgd = legend('e','a','s','t','x'); lgd.Location = 'southwest'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
subplot(1,2,2);
plot(errors');
title('Error var \sigma_e'); xlabel('Filter order'); ylabel('Error'); 
xlim([-2 length(filt_orders)]); 
lgd = legend('e','a','s','t','x'); lgd.Location = 'northwest'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
saveas(gcf,'CW4_q452_rp','epsc');






%% 4.5.3 Different sampling rate


% Adaptation
fs = 16000;
N = 1000;

% Resample audio
e_samp_l = resample(e, fs, info(1).SampleRate);
a_samp_l = resample(a, fs, info(1).SampleRate);

e_samp_l = e_samp_l(1:end);
e_samp_l= zscore(e_samp_l);
a_samp_l = a_samp_l(1:end);
a_samp_l= zscore(a_samp_l);


%% R_p etc.
Nw_max = 10;
filt_orders = 1:Nw_max;
lr = 0.001;

R_p = zeros(5, length(filt_orders));
errors = zeros(5, length(filt_orders));

w = 1;
letter = e_samp_l; lr = 0.01;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(1,w) = var(z_e);
    
    var_x = var(letter);
    var_e = var(z_e);
    R_p(1, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = a_samp; lr = 0.001;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(2,w) = var(z_e);
    
    var_x = var(letter);
    var_e = var(z_e);
    R_p(2, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = s_samp; lr = 0.005;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(3,w) = var(z_e);

    var_x = var(letter);
    var_e = var(z_e);
    R_p(3, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = t_samp; lr = 0.001;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(4,w) = var(z_e);

    var_x = var(letter);
    var_e = var(z_e);
    R_p(4, w) = 10*log10(var_x^2 / var_e^2);
    
    w= w+1;
end
w = 1;
letter = x_samp; lr = 0.001;
for filt_order = filt_orders
    [z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
    errors(5,w) = var(z_e);

    var_x = var(letter);
    var_e = var(z_e);
    R_p(5, w) = 10*log10(var_x / var_e);
    
    w= w+1;
end

% Normalize
for w=1:5
    R_p(w, :) = R_p(w, :)./max(R_p(w, :));
%     errors(i, :) = errors(i, :)./max(errors(i, :));
end

figure('PaperPosition', [0 0 30 10]);
subplot(1,2,1);
plot(R_p'); 
title('R_p by filter order'); xlabel('Filter order'); ylabel('R_p'); xlim([-2 length(filt_orders)]);
lgd = legend('e','a','s','t','x'); lgd.Location = 'southwest'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
subplot(1,2,2);
plot(errors');
title('Error var \sigma_e'); xlabel('Filter order'); ylabel('Error'); 
xlim([-2 length(filt_orders)]); 
lgd = legend('e','a','s','t','x'); lgd.Location = 'northwest'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
saveas(gcf,'CW4_q453_rp_N','epsc');



%% Quasi stationary signals
% Plot mean and stationary in window sizes for 'e' and 'a' at diff sampling
% freqs

window_size = 200;
mean_stat_e1 = zeros(1, length(window_size:window_size:length(e_samp)));
var_stat_e1 = zeros(1, length(window_size:window_size:length(e_samp)));
i=1;
for w = window_size:window_size:length(e_samp)
    mean_stat_e1(i) = mean(e_samp(w-window_size+1:w)); 
    var_stat_e1(i) = var(e_samp(w-window_size+1:w)); 
    i=i+1;
end
mean_stat_e2 = zeros(1, length(window_size:window_size:length(e_samp_l)));
var_stat_e2 = zeros(1, length(window_size:window_size:length(e_samp_l)));
i=1;
for w = window_size:window_size:length(e_samp_l)
    mean_stat_e2(i) = mean(e_samp_l(w-window_size+1:w)); 
    var_stat_e2(i) = var(e_samp_l(w-window_size+1:w)); 
    i=i+1;
end
mean_stat_a1 = zeros(1, length(window_size:window_size:length(a_samp)));
var_stat_a1 = zeros(1, length(window_size:window_size:length(a_samp)));
i=1;
for w = window_size:window_size:length(a_samp)
    mean_stat_a1(i) = mean(a_samp(w-window_size+1:w)); 
    var_stat_a1(i) = var(a_samp(w-window_size+1:w)); 
    i=i+1;
end
mean_stat_a2 = zeros(1, length(window_size:window_size:length(a_samp_l)));
var_stat_a2 = zeros(1, length(window_size:window_size:length(a_samp_l)));
i=1;
for w = window_size:window_size:length(a_samp_l)
    mean_stat_a2(i) = mean(a_samp_l(w-window_size+1:w)); 
    var_stat_a2(i) = var(a_samp_l(w-window_size+1:w)); 
    i=i+1;
end
figure('PaperPosition', [0 0 30 20]);   
subplot(2,2,1); hold on    % e1
plot(window_size:window_size:length(e_samp), mean_stat_e1); plot(window_size:window_size:length(e_samp),var_stat_e1);
title('Mean & var "e"'); xlabel('Index n'); ylabel('Stat magnitude');
lgd = legend('mean', 'var'); lgd.Location = 'northeast'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
subplot(2,2,2); hold on    % e2
plot(window_size:window_size:length(e_samp_l), mean_stat_e2); plot(window_size:window_size:length(e_samp_l),var_stat_e2);
title('Mean & var undersampled "e"'); xlabel('Index n'); ylabel('Stat magnitude');
lgd = legend('mean', 'var'); lgd.Location = 'northeast'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
subplot(2,2,3); hold on    % e1
plot(window_size:window_size:length(a_samp), mean_stat_a1); plot(window_size:window_size:length(a_samp),var_stat_a1);
title('Mean & var "a"'); xlabel('Index n'); ylabel('Stat magnitude');
lgd = legend('mean', 'var'); lgd.Location = 'northeast'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);
subplot(2,2,4); hold on    % a2
plot(window_size:window_size:length(a_samp_l), mean_stat_a2); plot(window_size:window_size:length(a_samp_l),var_stat_a2);
title('Mean & var undersampled "a"'); xlabel('Index n'); ylabel('Stat magnitude');
lgd = legend('mean', 'var'); lgd.Location = 'northeast'; lgd.FontSize = 15;
set(gca, 'Fontsize', 20);

saveas(gcf,'CW4_q453_stat','epsc');



%% 4.5.3 Stat signal 1000:2000

% Adaptation
fs = 16000;
i_start = 1000;
N = 1000;

% Resample audio
e_stat = resample(e, fs, info(1).SampleRate);
a_stat = resample(a, fs, info(1).SampleRate);

e_stat = e_stat(1:i_start+N);
e_stat= zscore(e_stat);
a_stat = a_stat(1:i_start+N);
a_stat= zscore(a_stat);


figure('PaperPosition', [0 0 25 8]);   
subplot(1,2,1); hold on
letter = a_stat; lr = 0.005;
filt_order = 2;
[z_hat, z_e, z_w] = lms_adap(letter(1:N), lr, filt_order, 0);
p=plot(1:N, z_w','color','g', 'LineWidth',2);
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[z_hat, z_e, z_w] = lms_adap(letter(i_start+1:i_start+N), lr, filt_order, 0);
p=plot(1:N, z_w', 'color','r');
set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
title('Weights undersampled "a"'); xlabel('Iteration N'); ylabel('Weight');
lgd = legend('1:1000','1001:2000'); lgd.Location = 'southeast';
set(gca, 'Fontsize', 20);
hold off

subplot(1,2,2); hold on
time_a = new_duration(2)*(i_start/N:1/N:(i_start+N-1)/N);
plot(time_a, a_stat(i_start+1:i_start+N)); plot(time_a, z_hat); 
title('Audio "a"'); xlabel('Time (s)'); ylabel('Audio Sig');
lgd = legend('Og a', 'LMS a'); lgd.FontSize = 12;
set(gca, 'Fontsize', 20);
hold off
saveas(gcf, 'CW4_q453_w', 'epsc');



