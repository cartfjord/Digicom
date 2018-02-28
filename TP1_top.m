close all

%Generate the PSS signals (pss0_t,pss1_t,pss2_t)
pss

% Do section 1.1 work here
figure(1);
subplot(3,1,1);
plot(real(pss0_t));
xlabel('sample n');
grid on;
title('I-Channel');

subplot(3,1,2);
plot(imag(pss0_t));
xlabel('sample n');
title('Q-Channel');
grid on;

subplot(3,1,3);
plot(abs(pss0_t));
xlabel('sample n');
title('Magnitude');
grid on;

% plot_handle = annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'PSS0 Sequence', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center');
% 
% plot_handle.FontSize = 12;
% plot_handle.FontWeight = 'bold';

figure(2);

% Do Sections 1.2 and 1.3 work hre
Fs = 15.36e6;
f_psd = -Fs/2:Fs/length(pss1_t):Fs/2;
f_psd = f_psd((end/2 - end/8) : (end/2 + end/8 - 1));

pss0_psd = 20*log10(abs(fftshift(fft(pss2_t))));
pss0_psd = pss0_psd((end/2 - end/8) : (end/2 + end/8 - 1));
plot(f_psd,pss0_psd);
title('PSS2 power spectrum');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
grid on;


figure(3);
subplot(3,1,1);
[xcorrel, lag] = xcorr(pss0_t);
plot(lag,abs(xcorrel));
title('autocorrelation pss0');

subplot(3,1,2);
[xcorrel, lag] = xcorr(pss1_t);
plot(lag,abs(xcorrel));
title('autocorrelation pss1');

subplot(3,1,3);
[xcorrel, lag] = xcorr(pss2_t);
plot(lag,abs(xcorrel));
title('autocorrelation pss2');


figure(4);
subplot(3,1,1);
[xcorrel01, lag] = xcorr(pss0_t,pss1_t);
plot(lag,abs(xcorrel01));
title('crosscorrelation between pss0 and pss1');

subplot(3,1,2);
[xcorrel02, lag] = xcorr(pss0_t,pss2_t);
plot(lag,abs(xcorrel02));
title('crosscorrelation between pss0 and pss2');

subplot(3,1,3);
[xcorrel12, lag] = xcorr(pss1_t,pss2_t);
plot(lag,abs(xcorrel12));
title('crosscorrelation between pss1 and pss2');

fprintf('ratio of signal energy pss0 and pss1: %.02f dB\n',-20*log10(abs(xcorrel01(1024))));
% 17.79 dB
fprintf('ratio of signal energy pss0 and pss2: %.02f dB\n',-20*log10(abs(xcorrel02(1024))));
% 8.30 dB
fprintf('ratio of signal energy pss0 and pss2: %.02f dB\n',-20*log10(abs(xcorrel12(1024))));
% 17.79 dB

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

clear s;

figure(5);
subplot(2,3,1);
%plot an approximation to the power spectrum
f = linspace(-7.68e6,7.68e6,153600*8);
plot(f,20*log10(abs(fftshift(fft(s_796)))))
axis([-7.68e6 7.68e6 100 150])
title('Power spectrum of f_c = 796 MHz signal');
xlabel('f');
ylabel('Magnitude [dB]');

subplot(2,3,4);
t = linspace(0,80e-3,length(s_796));
plot(t,mag2db(abs(s_796)));
%bitsra(length(s_796),3);
title('Time plot of f_c = 796 MHz signal');
xlabel('t');
ylabel('Magnitude [dB]');

subplot(2,3,2);
%plot an approximation to the power spectrum
f = linspace(-7.68e6,7.68e6,153600*8);
plot(f,20*log10(abs(fftshift(fft(s_806)))))
axis([-7.68e6 7.68e6 100 150])
title('Power spectrum of f_c = 806 MHz signal');
xlabel('f');
ylabel('Magnitude [dB]');

subplot(2,3,5);
t = linspace(0,80e-3,length(s_796));
plot(t,mag2db(abs(s_806)));
title('Time plot of f_c = 806 MHz signal');
xlabel('t');
ylabel('Magnitude [dB]');

subplot(2,3,3);
%plot an approximation to the power spectrum
f = linspace(-7.68e6,7.68e6,153600*8);
plot(f,20*log10(abs(fftshift(fft(s_816)))))
axis([-7.68e6 7.68e6 100 150])
title('Power spectrum of f_c = 816 MHz signal');
xlabel('f');
ylabel('Magnitude [dB]');

subplot(2,3,6);
t = linspace(0,80e-3,length(s_816));
plot(t,mag2db(abs(s_816)));
title('Time plot of f_c = 816 MHz signal');
xlabel('t');
ylabel('Magnitude [dB]');
% 
% 
% [i, Nf] = max_pss(s_796);
% fprintf('pss: %d, Nf = %d\n', i, Nf);
% switch i
%     case 0
%         freq_offset_est(s_796, pss0_t, Nf);
%     case 1
%         freq_offset_est(s_796, pss1_t, Nf);
%     case 2
%         freq_offset_est(s_796, pss2_t, Nf);
% end
% 
% 
[i, Nf] = max_pss(s_806);
fprintf('pss: %d, Nf = %d\n', i, Nf);
switch i
    case 0
        freq_offset_est(s_806, pss0_t, Nf);
    case 1
        freq_offset_est(s_806, pss1_t, Nf);
    case 2
        freq_offset_est(s_806, pss2_t, Nf);
end
% 
% 
% [i, Nf] = max_pss(s_816);
% fprintf('pss: %d, Nf = %d\n', i, Nf);
% switch i
%     case 0
%         freq_offset_est(s_816, pss0_t, Nf);
%     case 1
%         freq_offset_est(s_816, pss1_t, Nf);
%     case 2
%         freq_offset_est(s_816, pss2_t, Nf);
% end
% 
% %Find the frequency offset
% 
% 
