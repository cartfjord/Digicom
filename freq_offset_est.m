function [ f_offset ] = freq_offset_est(signal, pss_i, Nf)
    %%Frequency offset estimator
    
    %To make sure the signals are generated before running this func. 
    pss; 
    
    figure;
    %suptitle('Frequncy offset'); %Requires Biometrics toolbox in Matlab.
    DELTA_F = 100;
    

    Fs = 15.36e6;
    f_min = -7500;
    f_max = 7500;
    
    m = f_min:DELTA_F:f_max;
    Y = zeros(1,length(m)); 
  
    t = 0:(1/Fs):(1023/Fs); 
    
    %figure;
    for sample_offset = 0:76800:76800*15 % 5 ms steps
        for j = 1:length(m)
            Y(j) = Y(j) + abs(sum(exp(-2*pi*1i*m(j).*t).*conj(pss_i).*signal((Nf + sample_offset):(Nf + sample_offset + 1023)).')).^2;
        end
    end

    plot(m,Y);
    
    title('Most likely frequency offset');
    xlabel('Frequency offset');
    ylabel('Likelihood');
    
    grid on;
    [~, f_offset] = max(Y);
    fprintf('Detected offset: %d Hz\n',m(f_offset));
end

