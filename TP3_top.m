cell_id_796 = 288;
cell_id_806 = 365;
cell_id_816 = 421;

Ncp_796 = 0;
Ncp_806 = 0;
Ncp_816 = 0;

N_RB_DL = 50;

% Read in sample file (example here) 
fd = fopen('signal806_2.dat','r') ;
s = fread(fd,153600*2*8,'int16') ; 
fclose(fd) ; 
s_806 = s(1:2:end) + sqrt(-1)*s(2:2:end) ; 

if(~exist('lte_gold_table_796','var') || ...
   ~exist('lte_gold_table_806','var') || ...
   ~exist('lte_gold_table_816','var') )

    lte_gold_table_796 = lte_rs_gold(Ncp_796, cell_id_796);
    lte_gold_table_806 = lte_rs_gold(Ncp_806, cell_id_806);
    lte_gold_table_816 = lte_rs_gold(Ncp_816, cell_id_816);

end

l = 0; % Symbol within a slot
Ns = 1; % Slot index
p = 0   ; % Antenna

dl_cs_rs_796 = lte_rs(Ncp_796, cell_id_796, l, Ns, N_RB_DL, lte_gold_table_796);
dl_cs_rs_806 = lte_rs(Ncp_806, cell_id_806, l, Ns, N_RB_DL, lte_gold_table_806);
dl_cs_rs_816 = lte_rs(Ncp_816, cell_id_816, l, Ns, N_RB_DL, lte_gold_table_816);


%% Generation of ks
nu = 3*(l ~= p);

nu_shift_796 = mod(cell_id_796, 6);
nu_shift_806 = mod(cell_id_806, 6);
nu_shift_816 = mod(cell_id_816, 6);

n = 0 : (2 * N_RB_DL - 1);

k_796 = 6 * n + mod((nu + nu_shift_796), 6);
k_806 = 6 * n + mod((nu + nu_shift_806), 6);
k_816 = 6 * n + mod((nu + nu_shift_816), 6);

%% continue

Nf_0 = 72406; %38048;
Nf_1 = 85713; %22494; %8913; 
Nf_2 = 47954; %33339;


%pbch_796 = fft(s_796(Nf_0 + 1096 : Nf_0 + 1096 + 1023));
pbch_806_1 = fftshift(fft(s_806(Nf_1 + 1*(1024 + 72) + 1*76800 + (1:1024)))); % 163610 - 164633
pbch_806_2 = fftshift(fft(s_806(Nf_1 + 2*(1024 + 72) + 1*76800 + (1:1024)))); % 164706 - 165729
pbch_806_3 = fftshift(fft(s_806(Nf_1 + 3*(1024 + 72) + 1*76800 + (1:1024)))); % 165802 - 166825
pbch_806_4 = fftshift(fft(s_806(Nf_1 + 4*(1024 + 72) + 1*76800 + (1:1024)))); % 166898 - 167921
after_pbch = fftshift(fft(s_806(Nf_1 + 5*(1024 + 72) + 1*76800 + (1:1024)))); % 167994 - 169017

pbch_806_1_600 = pbch_806_1(512 - 299 : 512 + 300);
pbch_806_2_600 = pbch_806_2(512 - 299 : 512 + 300);
pbch_806_3_600 = pbch_806_3(512 - 299 : 512 + 300);
pbch_806_4_600 = pbch_806_4(512 - 299 : 512 + 300);
after_pbch_600 = after_pbch(512 - 299 : 512 + 300);

%%Extract channel estimates
chan_est_symb_0 = transpose(pbch_806_1(k_806)).*dl_cs_rs_806;
chan_est_symb_1 = transpose(pbch_806_2(k_806)).*dl_cs_rs_806;
chan_est_symb_2 = transpose(pbch_806_2(k_806)).*dl_cs_rs_806;
chan_est_symb_3 = transpose(pbch_806_2(k_806)).*dl_cs_rs_806;
chan_est_symb_after = transpose(after_pbch(k_806)).*dl_cs_rs_806;

%%Interpolate channel estimates
chan_est_symb_0_int = interp1(k_806, chan_est_symb_0, 1:600, 'linear','extrap');
chan_est_symb_1_int = interp1(k_806, chan_est_symb_1, 1:600, 'linear','extrap');
chan_est_symb_2_int = interp1(k_806, chan_est_symb_2, 1:600, 'linear','extrap');
chan_est_symb_3_int = interp1(k_806, chan_est_symb_3, 1:600, 'linear','extrap');
chan_est_symb_after_int = interp1(k_806, chan_est_symb_after, 1:600, 'linear','extrap');



figure;

%Plot channel estimate over time

surf([  abs(chan_est_symb_0_int);...
        abs(chan_est_symb_1_int);...
        abs(chan_est_symb_2_int);...
        abs(chan_est_symb_3_int);...
        abs(chan_est_symb_after_int)
        ]);
xlabel('Interpolated subcarrier index');
ylabel('Time [symbol index]');
zlabel('|H(f)|');
title('Channel estimate for PBCH over time');

surf([  abs(ifft(chan_est_symb_0_int));...
        abs(ifft(chan_est_symb_1_int));...
        abs(ifft(chan_est_symb_2_int));...
        abs(ifft(chan_est_symb_3_int));...
        abs(ifft(chan_est_symb_after_int))
        ]);
%stem(abs(ifft(chan_est_symb_3_int)));
xlim([1 300]);
grid on;
xlabel('Sample [n]');
ylabel('Time [symbol index]');
zlabel('|h(t)|');
title('Channel impulse response');


