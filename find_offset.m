function Nf = find_offset(signal)
    
    window_len = floor((size(signal) + 1) / 8);
    window_len = window_len(1);

    signal = signal';
    signal = [signal 0];
    added_peaks = zeros(1,window_len);
    
    for i = 1:window_len:7*window_len
        to_add = signal(i:(i+window_len-1));
        added_peaks = added_peaks + to_add;
    end
    
    Nf = abs(added_peaks);
    
end
