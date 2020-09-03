% ASP CW 4: 4.6
clear all;
close all;
rng(0);

%% 4.6

% AR COEFFS
N = 1000;
wgn = randn(1, N);
a = [1, 0.9, 0.2];
b = [1];
x = filter(b, a, wgn);

% MSE in dB
letter = x; lr = 0.01; filt_order = length(a)-1; is_gear_shift = 1;
% sign error
[z_hat, z_e_e, z_w_e] = sign_lms(letter, lr, filt_order, 1,0, is_gear_shift);
% sign reg
[z_hat, z_e_r, z_w_r] = sign_lms(letter, lr, filt_order, 0,1, is_gear_shift);
% sign sign
[z_hat, z_e, z_w] = sign_lms(letter, lr, filt_order, 1,1, is_gear_shift);
% LMS
[z_hat, z_e_lms, z_w_lms] = lms_adap(letter, lr, filt_order, is_gear_shift);

% Normalize 
z_e_e = 10*log10(abs(z_e_e)./max(z_e_e));
z_e_r = 10*log10(abs(z_e_r)./max(z_e_r));
z_e = 10*log10(abs(z_e)./max(z_e));
z_e_lms = 10*log10(abs(z_e_lms)./max(z_e_lms));

z_w_e = (abs(z_w_e))'; %./max(z_e_e)));
z_w_r = (abs(z_w_r))'; %./max(z_e_r)));
z_w = (abs(z_w))'; %./max(z_e)));
z_w_lms = (abs(z_w_lms))'; %./max(z_e_lms)));


c = 'r';
figure('PaperPosition', [0 0 45 8]);   
subplot(1,4,1); plot((z_e_e), 'color',c); 
title('Sign error');
set(gca, 'Fontsize', 20); ylabel('MSE (dB)');
subplot(1,4,2); plot((z_e_r), 'color',c);
title('Sign regressor');
set(gca, 'Fontsize', 20); xlabel('Iteration n');
subplot(1,4,3); plot((z_e), 'color',c); 
title('Sign sign');
set(gca, 'Fontsize', 20);
subplot(1,4,4); plot((z_e_lms), 'color',c); 
title('Classic LMS');
set(gca, 'Fontsize', 20);
sgtitle('Sign LMS Learning Curves AR coeffs', 'Fontsize', 20);

saveas(gcf, 'CW4_q461_ar','epsc');


c = 'r';
figure('PaperPosition', [0 0 40 7]);   
subplot(1,4,1); plot((z_w_e), 'color',c);
title('Sign error');
set(gca, 'Fontsize', 20); ylabel('MSE (dB)');
txt = strcat('w(1)=', num2str(z_w_e(end,1)));
txt2 = strcat('w(2)=', num2str(z_w_e(end,2)));
text(N/3, max(z_w_e(end,:))/2,txt,'FontSize',15);
text(N/3, max(z_w_e(end,:))/2-0.2, txt2,'FontSize',15);
subplot(1,4,2); plot((z_w_r), 'color',c);
title('Sign regressor');
set(gca, 'Fontsize', 20); xlabel('Iteration n');
txt = strcat('w(1)=', num2str(z_w_r(end,1)));
txt2 = strcat('w(2)=', num2str(z_w_r(end,2)));
text(N/3, max(z_w_r(end,:))/2,txt,'FontSize',15);
text(N/3, max(z_w_r(end,:))/2-0.2, txt2,'FontSize',15);
subplot(1,4,3); plot((z_w), 'color',c); 
title('Sign sign');
set(gca, 'Fontsize', 20);
txt = strcat('w(1)=', num2str(z_w(end,1)));
txt2 = strcat('w(2)=', num2str(z_w(end,2)));
text(N/3, max(z_w(end,:))/2,txt,'FontSize',15);
text(N/3, max(z_w(end,:))/2-0.2, txt2,'FontSize',15);
subplot(1,4,4); plot((z_w_lms), 'color',c); 
title('Classic LMS');
set(gca, 'Fontsize', 20);
txt = strcat('w(1)=', num2str(z_w_lms(end,1)));
txt2 = strcat('w(2)=', num2str(z_w_lms(end,2)));
text(N/3, max(z_w_lms(end,:))/2,txt,'FontSize',15);
text(N/3, max(z_w_lms(end,:))/2-0.2, txt2,'FontSize',15);
sgtitle(['Weights Sign LMS AR coeffs, lr=',num2str(lr)], 'Fontsize', 20);

saveas(gcf, 'CW4_q461_arw','epsc');


% R_p
R_p = zeros(1,4);
% errors = zeros(5, length(filt_orders));
% 
% errors(1) = var(z_e);
var_x = var(letter);
var_e = var(z_e_e);
R_p(1) = 10*log10(var_x^2 / var_e^2);
var_e = var(z_e_r);
R_p(2) = 10*log10(var_x^2 / var_e^2);
var_e = var(z_e);
R_p(3) = 10*log10(var_x^2 / var_e^2);
var_e = var(z_e_lms);
R_p(4) = 10*log10(var_x^2 / var_e^2);

figure('PaperPosition', [0 0 20 7]);   
bar(R_p); 
title('R_p Sign LMS AR'); xlabel('Algorithm'); ylabel('R_p (dB)');
set(gca, 'Fontsize', 20);
set(gca,'XTickLabel',{'Error','Regressor','Sign','LMS'});
saveas(gcf, 'CW4_q461_arrp','epsc');



%% SPEECH 
% Audio: come in 2 identical columns, drop 1
e = audioread('e.wav');
info = audioinfo('e.wav');
% Delete silence
window_size = 300;
e = del_silence(e, window_size, 1);
e = zscore(e);

% MSE in dB
letter = e(1000:2000); lr = 0.001; filt_order = 2; is_gear_shift=0;
% sign error
[z_hat, z_e_e, z_w] = sign_lms(letter, lr, filt_order, 1,0, is_gear_shift);
% sign reg
[z_hat, z_e_r, z_w] = sign_lms(letter, lr, filt_order, 0,1, is_gear_shift);
% sign sign
[z_hat, z_e, z_w] = sign_lms(letter, lr, filt_order, 1,1, is_gear_shift);
% LMS
[z_hat, z_e_lms, z_w] = lms_adap(letter, lr, filt_order, is_gear_shift);

% Normalize 
z_e_e = 10*log10(abs(z_e_e)); %./max(z_e_e)));
z_e_r = 10*log10(abs(z_e_r)); %./max(z_e_r)));
z_e = 10*log10(abs(z_e)); %./max(z_e)));
z_e_lms = 10*log10(abs(z_e_lms)); %./max(z_e_lms)));

figure('PaperPosition', [0 0 45 10]);   
subplot(1,4,1); plot((z_e_e)); 
title('Sign error');
set(gca, 'Fontsize', 20); ylabel('MSE (dB)');
subplot(1,4,2); plot((z_e_r));
title('Sign regressor');
set(gca, 'Fontsize', 20); xlabel('Iteration n');
subplot(1,4,3); plot((z_e)); 
title('Sign sign');
set(gca, 'Fontsize', 20);
subplot(1,4,4); plot((z_e_lms)); 
title('Classic LMS');
set(gca, 'Fontsize', 20);
sgtitle('Learning Curves Sign LMS "e"', 'Fontsize', 20);

saveas(gcf, 'CW4_q461','epsc');


% R_p
R_p = zeros(1,4);
var_x = var(letter);
var_e = var(z_e_e);
R_p(1) = 10*log10(var_x^2 / var_e^2);
var_e = var(z_e_r);
R_p(2) = 10*log10(var_x^2 / var_e^2);
var_e = var(z_e);
R_p(3) = 10*log10(var_x^2 / var_e^2);
var_e = var(z_e_lms);
R_p(4) = 10*log10(var_x^2 / var_e^2);

figure('PaperPosition', [0 0 20 7]);   
bar(R_p); 
title('R_p Sign LMS "a"'); xlabel('Algorithm'); ylabel('R_p (dB)');
set(gca, 'Fontsize', 20);
set(gca,'XTickLabel',{'Error','Regressor','Sign','LMS'});
saveas(gcf, 'CW4_q461_a_rp','epsc');
