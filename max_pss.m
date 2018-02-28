function [i, Nf] = max_pss(signal)
    pss; % Make sure the signals have been generated.
    
    % Implementation of filters
    xcorr0 = conv(signal,fliplr(conj(pss0_t)));
    xcorr1 = conv(signal,fliplr(conj(pss1_t)));
    xcorr2 = conv(signal,fliplr(conj(pss2_t)));
    
    % Find the PSS with the highest peaks
    [~,i] = max([max(xcorr0), max(xcorr1), max(xcorr2)]);
    i = i - 1;
    
    %Find the highest peak
    switch i
        case 0
            [~, Nf] = max(xcorr0);
        case 1
            [~, Nf] = max(xcorr1);
        case 2
            [~, Nf] = max(xcorr2);
    end
    
    %Since we know the higest peak and Fs,
    %we can find the first peak:
    Nf = Nf - 1023;
    
    while(Nf > 76800) 
        Nf = Nf - 76800;
    end
    
    %Plots
    figure
    
    %suptitle(['PSS = ' num2str(i), ' Nf = ', num2str(Nf)]);
    subplot(3,1,1);
    plot(abs(xcorr0));
    title('Matched filter applied to PSS_0');
    xlabel('Sample [n]');
    grid on;
    subplot(3,1,2);
    plot(abs(xcorr1));
    title('Matched filter applied to PSS_1');
    xlabel('Sample [n]');
    grid on;
    subplot(3,1,3);
    plot(abs(xcorr2));
    title('Matched filter applied to PSS_2');
    xlabel('Sample [n]');
    grid on;
end
