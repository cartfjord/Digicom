function [dl_cs_rs] = lte_rs(Ncp,Nid_cell,l,Ns,N_RB_DL,lte_gold_table)

% includes conjugate transpose for channel estimation at RX
qpsk(1) = (1-sqrt(-1))/sqrt(2);
qpsk(2) = (-1-sqrt(-1))/sqrt(2);
qpsk(3) = (1+sqrt(-1))/sqrt(2);
qpsk(4) = (-1+sqrt(-1))/sqrt(2);

mprime = 110 - N_RB_DL;
k=0;

for m=0 : (2*N_RB_DL - 1),

    mprime_dword     = uint32(floor(mprime/16));
    mprime_qpsk_symb = uint32(rem(mprime,16));

    % this is r_mprime from 3GPP 36-211 6.10.1.2
    dl_cs_rs(1+k) = qpsk(1+rem(lte_gold_table(1+Ns,1+l,1+mprime_dword)/(2*mprime_qpsk_symb),4));
    fprintf('Ns %d, l %d, m %d,mprime_dword %d, mprime_qpsk_symbol %d\n',Ns,l,m,mprime_dword,mprime_qpsk_symb);
    fprintf('index = %d (k %d)\n',rem(lte_gold_table(1+Ns,1+l,1+mprime_dword)/(2*mprime_qpsk_symb),4),k);

    mprime=mprime+1;

    if (m<4)
      fprintf('Ns %d l %d output[%d] = (%f,%f)\n',Ns,l,k,real(dl_cs_rs(1+k)),imag(dl_cs_rs(1+k)));
    end
    
    k=k+1;
    
end