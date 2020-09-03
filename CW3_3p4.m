% ASP CW3 3.4
clear all
close all

%% 3.4.1 

random_phone = [0 2 0 floor(rand(1, 8)*10)];
nums = random_phone(4:5);
N = length(nums); % length(random_phone);

k = {0,1,2,3,4,5,6,7,8,9};
v = {[1336 941], [1209 697], [1336 697], [1477 697], [1209 770], [1336 770], [1477 770], [1209 852], [1336 852], [1477 852]};
freq_map = containers.Map(k, v);

f = zeros((2*N-1),2);
fs = 32768;
% t1 = 0.25*fs;

% Make dial tone 1 0 1 0 ... 1
x1 = (1/fs):(1/fs):0.25;
y=[];
pause = zeros(1, length(x1));         
for num=nums
    f = freq_map(num);
    y_temp = sin(2*pi*f(1)*x1) + sin(2*pi*f(2)*x1);
    if num==nums(end)
        y=[y y_temp];
    else
        y=[y y_temp pause];
    end
end
% f = freq_map(random_phone(num+1));
% y_temp = sin(2*pi*f(1)*x1) + sin(2*pi*f(2)*x1);
% y=[y y_temp];

t = (1/fs):(1/fs):(N+1)*0.25;
figure('PaperPosition', [0 0 25 8]); 
plot(t,y);
title(['Dial tones ', num2str(nums)]);
xlabel('time (s)'); ylabel('Amplitude');
set(gca, 'Fontsize', 18);
saveas(gcf, 'CW3_q341', 'epsc');


%% Nyquist rate: find fundamental freq

gcd_dials = zeros(1, length(k));
for num=0:9
    f = freq_map(num);
    gcd_dials(num+1) = gcd(f(1), f(2));
end
f0 = max(max(gcd_dials) * freq_map(find(gcd_dials(gcd_dials==max(gcd_dials)))-1));





%% 3.4.2
% get full y
nums = random_phone;
y=[];
pause = zeros(1, length(x1));         
for num=nums
    f = freq_map(num);
    y_temp = sin(2*pi*f(1)*x1) + sin(2*pi*f(2)*x1);
    if num==nums(end)
        y=[y y_temp];
    else
        y=[y y_temp pause];
    end
end

figure;
n_sect = length(x1);
n_overl = 0;
spectrogram(y, hanning(n_sect), n_overl, n_sect, fs);
title(['Dial Tone Spectrogram, ', num2str(random_phone)]);
xlim([0 N]);
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW3_q342', 'epsc');



%% 3.4.4 REPEAT LOL

y3 = y + randn(size(y));
y2 = y + 0.5*randn(size(y));
y1 = y + 0.1*randn(size(y));

% 1
N_full = length(y2);
t = linspace(0, 5.25, N_full);

figure('PaperPosition', [0 0 30 8]); 
subplot(1,3,1); plot(t,y1);
title(['Min. Noisy Dial']);
ylabel('Amplitude');
set(gca, 'Fontsize', 20);
subplot(1,3,2); plot(t, y2);
title('Med. Noisy Dial');
xlabel('time (s)'); 
set(gca, 'Fontsize', 20);
subplot(1,3,3); plot(t, y3);
title('Max. Noisy Dial');
set(gca, 'Fontsize', 20);
saveas(gcf, 'CW3_q345_1', 'epsc');

% 2
n_sect = length(x1);
n_overl = 0;
figure('PaperPosition', [0 0 30 8]); 
subplot(1,3,1); spectrogram(y1, hanning(n_sect), n_overl, n_sect, fs);
title(['Min. Noisy Spectrogram']);
xlim([0 N]); 
set(gca, 'Fontsize', 18);
subplot(1,3,2); spectrogram(y2, hanning(n_sect), n_overl, n_sect, fs);
title(['Med. Noisy Spectr.']);
xlim([0 N]); ylabel('');
set(gca, 'Fontsize', 18);
subplot(1,3,3); spectrogram(y3, hanning(n_sect), n_overl, n_sect, fs);
title(['Max. Noisy Spectr.']);
xlim([0 N]); ylabel('');
set(gca, 'Fontsize', 18);
saveas(gcf, 'CW3_q345_2', 'epsc');


