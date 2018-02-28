close all;
clear all;

pss;

Nf_0 = 72406; %38048;
Nf_1 = 22495; %22495; %8913;
Nf_2 = 47954; %33339;

offset_0 = 834; %400
offset_1 = 790; %200
offset_2 = 944; %200

cell_0 = 2;
cell_1 = 2;
cell_2 = 1;

% Read in sample file (example here)
%fd = fopen('signal796.dat','r') ; 
fd = fopen('center_796_2.dat','r') ;
s = fread(fd,153600*2*8,'int16') ; 
fclose(fd) ; 
s_796 = s(1:2:end) + sqrt(-1)*s(2:2:end) ; 

% Read in sample file (example here)
%fd = fopen('signal806.dat','r') ; 
fd = fopen('center_806_2.dat','r') ;
s = fread(fd,153600*2*8,'int16') ; 
fclose(fd) ; 
s_806 = s(1:2:end) + sqrt(-1)*s(2:2:end) ; 

% Read in sample file (example here)
%fd = fopen('signal816.dat','r') ; 
fd = fopen('center_816_2.dat','r') ;
s = fread(fd,153600*2*8,'int16') ; 
fclose(fd) ; 
s_816 = s(1:2:end) + sqrt(-1)*s(2:2:end) ; 

Fs = 15.36e6;
fft_len = 1024;
f = -Fs/2:Fs/fft_len:Fs/2 -1;

t = 0:(1/Fs):(1023/Fs);

apply_freq_offset = 1;

if(apply_freq_offset)
    s_796 = transpose(s_796).*exp(-2*pi*1i*offset_0.*(1:length(s_796))/Fs);
    s_806 = transpose(s_806).*exp(-2*pi*1i*offset_1.*(1:length(s_806))/Fs);
    s_816 = transpose(s_816).*exp(-2*pi*1i*offset_2.*(1:length(s_816))/Fs);
end

%figure;
%%subplot(1,2,1);

% pss_0_subframe0 = s_796(Nf_0 + (1:1024));
% sss_0_subframe0 = s_796((Nf_0 - 1024 - 72):(Nf_0 - 1024 - 72 + 1023));
% pss_0_subframe5 = s_796(Nf_0  + 76800 + (1:1024));
% sss_0_subframe5 = s_796((Nf_0 + 76800 - 1024 - 72):(Nf_0 + 76800 - 1024 - 72 + 1023));
% 
% PSS_0_subframe0 = fft(pss_0_subframe0, fft_len);
% SSS_0_subframe0 = fft(sss_0_subframe0, fft_len);
% PSS_0_subframe5 = fft(pss_0_subframe5, fft_len);
% SSS_0_subframe5 = fft(sss_0_subframe5, fft_len);
% 
% chan_0_response = PSS_0_subframe0 .* conj(pss0_f);
% 
% PSS_0_extract = [PSS_0_subframe0((1024-30):1024) PSS_0_subframe0(2:32)];
% 
% 
% stem(f,fftshift(abs(chan_0_response)));
% xlim([-5.5e5 5.5e5]);

% subplot(1,2,2);
% SSS_0_62 = [SSS_0_subframe0((512-30):512),SSS_0_subframe0(514:(514+30))];
% 
% 
% scatter(real(SSS_0_62),imag(SSS_0_62));
% axis('square');
% 
% figure;
% subplot(1,3,1);
% 
pss_1_subframe0 = s_806(Nf_1 + (1:1024));
sss_1_subframe0 = s_806((Nf_1 - 1024 - 72):(Nf_1 - 1024 - 72 + 1023));
pss_1_subframe5 = s_806(Nf_1  + 76800 + (1:1024));
sss_1_subframe5 = s_806((Nf_1 + 76800 - 1024 - 72):(Nf_1 + 76800 - 1024 - 72 + 1023));

PSS_1_subframe0 = fft(pss_1_subframe0, fft_len);
SSS_1_subframe0 = fft(sss_1_subframe0, fft_len);
PSS_1_subframe5 = fft(pss_1_subframe5, fft_len);
SSS_1_subframe5 = fft(sss_1_subframe5, fft_len);

figure;
subplot(1,2,1);
stem(f,fftshift(abs(SSS_1_subframe0)));
xlim([-5.5e5 5.5e5]);
grid on;
xlabel('Frequency [Hz]');
ylabel('|SSS|');
title('Six PRBSs around 0 Hz at PSS subframe 0, f_c = 806 MHz');
subplot(1,2,2);
plot(f,fftshift(abs(SSS_1_subframe0)));
grid on;
xlabel('Frequency [Hz]');
ylabel('|SSS|');
title('Entire OFDM symbol, f_c = 806 MHz');
figure;

chan_1_response = PSS_1_subframe0 .* conj(pss2_f);
chan_1_response_subframe5 = PSS_1_subframe5 .* conj(pss2_f);
PSS_1_extract = [PSS_1_subframe0((1024-30):1024) PSS_1_subframe0(2:32)];
subplot(2,1,1);
stem(f,fftshift(abs(chan_1_response)));
xlim([-5.5e5 5.5e5]);
grid on;
xlabel('Frequency [Hz]');
ylabel('|H_{pss}|');
title('Channel frequnecy response f_c = 806 MHz');

subplot(2,1,2);
chan_1_72_impulse_response = ifft([chan_1_response(1 : 36), chan_1_response(end - 36 : end)]);
chan_1_72_impulse_response = chan_1_72_impulse_response(1 : end/2);
stem(0:35, abs(chan_1_72_impulse_response));

grid on;
xlabel('n');
ylabel('|h_{pss}|');
title('Channel impulse response f_c = 806 MHz');

% subplot(1,3,2);
figure;
SSS_1_62_subframe0 = [SSS_1_subframe0((1024-30):1024) SSS_1_subframe0(2:32)];
SSS_1_62_subframe5 = [SSS_1_subframe5((1024-30):1024) SSS_1_subframe5(2:32)];

hold on;
scatter(real(SSS_1_62_subframe0),imag(SSS_1_62_subframe0));
scatter(real(SSS_1_62_subframe5),imag(SSS_1_62_subframe5));
hold off;

title('f_c = 806 MHz, uncompensated for channel');
axis('square');
xlim([-3e5, 3e5]);
ylim([-3e5, 3e5]);
grid on;
xlabel('I');
ylabel('Q');
legend('SSS subframe 0', 'SSS subframe 5');
% 
% subplot(1,3,3);

figure;
chan_1_62_subframe0 = [chan_1_response((1024-30):1024) chan_1_response(2:32)];
chan_1_62_subframe5 = [chan_1_response_subframe5(1024 - 30 : 1024) chan_1_response_subframe5(2 : 32)];

% 
SSS_1_62_subframe0_comp = SSS_1_62_subframe0.*conj(chan_1_62_subframe0);
SSS_1_62_subframe5_comp = SSS_1_62_subframe5.*conj(chan_1_62_subframe5);


hold on;
scatter(real(SSS_1_62_subframe0_comp),imag(SSS_1_62_subframe0_comp));
scatter(real(SSS_1_62_subframe5_comp),imag(SSS_1_62_subframe5_comp));
hold off;

title('f_c = 806 MHz, channel compensated');
axis('square');
xlim([-2e15, 2e15]);
ylim([-2e15, 2e15]);
grid on;
xlabel('I');
ylabel('Q');
legend('SSS subframe 0', 'SSS subframe 5');

figure;
subplot(1,3,1);

pss_2_subframe0 = s_816(Nf_2 + (1:1024));
sss_2_subframe0 = s_816((Nf_2 - 1024 - 72):(Nf_2 - 1024 - 72 + 1023));
pss_2_subframe5 = s_816(Nf_2  + 76800 + (1:1024));
sss_2_subframe5 = s_816((Nf_2 + 76800 - 1024 - 72):(Nf_2 + 76800 - 1024 - 72 + 1023));

PSS_2_subframe0 = fft(pss_2_subframe0, fft_len);
SSS_2_subframe0 = fft(sss_2_subframe0, fft_len);
PSS_2_subframe5 = fft(pss_2_subframe5, fft_len);
SSS_2_subframe5 = fft(sss_2_subframe5, fft_len);

chan_2_response_subframe0 = PSS_2_subframe0 .* conj(pss1_f);
chan_2_response_subframe5 = PSS_2_subframe5 .* conj(pss1_f);

PSS_2_extract = [PSS_2_subframe0((1024-30):1024) PSS_2_subframe0(2:32)];

stem(f,fftshift(abs(chan_2_response_subframe0)));
xlim([-5.5e5 5.5e5]);


subplot(1,3,2);
SSS_2_62_subframe0 = [SSS_2_subframe0((1024-30):1024) SSS_2_subframe0(2:32)];
SSS_2_62_subframe5 = [SSS_2_subframe5((1024-30):1024) SSS_2_subframe5(2:32)]; 

chan_2_62_subframe0 = [chan_2_response_subframe0((1024-30):1024) chan_2_response_subframe0(2:32)];
chan_2_62_subframe5 = [chan_2_response_subframe5((1024-30):1024) chan_2_response_subframe5(2:32)];

SSS_2_62_subframe0_comp = SSS_2_62_subframe0.*conj(chan_2_62_subframe0);
SSS_2_62_subframe5_comp = SSS_2_62_subframe5.*conj(chan_2_62_subframe5);

hold on;
scatter(real(SSS_2_62_subframe0),imag(SSS_2_62_subframe0));
scatter(real(SSS_2_62_subframe5),imag(SSS_2_62_subframe5)); %HACK
hold off;

title('816 MHz, uncompensated');
xlim([-1e5, 1e5]);
ylim([-1e5, 1e5]);
axis('square');

subplot(1,3,3);


hold on;

scatter(real(SSS_2_62_subframe0_comp),imag(SSS_2_62_subframe0_comp));
scatter(real(SSS_2_62_subframe5_comp),imag(SSS_2_62_subframe5_comp)); %HACK

hold off;

title('816 MHz, compensated');
axis('square');
xlim([-8e13, 8e13]);
ylim([-8e13, 8e13]);


% figure;
% subplot(1,3,1);
% plot(abs(ifft(chan_0_response)));
% title('796 MHz |H(f)|');
% subplot(1,3,2)
% plot(abs(ifft(chan_1_response)));
% title('806 MHz |H(f)|');
% subplot(1,3,3);
% plot(abs(ifft(chan_2_response_subframe0)));
% title('816 MHz |H(f)|');

sss; %Generate d0 and d5

%3. Data detection

SSS_1_62_subframe0_comp_sliced = 2*(real(SSS_1_62_subframe0_comp) > 0) - 1;
SSS_1_62_subframe5_comp_sliced = 2*(real(SSS_1_62_subframe5_comp) > 0) - 1;

SSS_2_62_subframe0_comp_sliced = 2*(real(SSS_2_62_subframe0_comp) > 0) - 1;
SSS_2_62_subframe5_comp_sliced = 2*(real(SSS_2_62_subframe5_comp) > 0) - 1;

%% Determine the most likely Cell ID
[cell_id_806_subframe0, cell_id_806_max_subframe0] = find_best_sss(SSS_1_62_subframe0_comp_sliced);
[cell_id_806_subframe5, cell_id_806_max_subframe5] = find_best_sss(SSS_1_62_subframe5_comp_sliced);

if(cell_id_806_max_subframe0 > cell_id_806_max_subframe5)
    cell_id_806 = cell_id_806_subframe0;
else
    cell_id_806 = cell_id_806_subframe5;
end

disp(cell_id_806);
