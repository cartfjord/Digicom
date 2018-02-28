function [sss_idx, sss_max] = find_best_sss(subframe)
sss;
% Exhausitve search can be optimizied, but performance is not really an
% issue when we are writing in Matlab in the first place.
sss_idx = 0;
sss_max = 0;

for i = 1:504   
   sss_max_tmp = sum(xcorr(d0(i,:),subframe));
   if(abs(sss_max_tmp) >= sss_max)
       sss_max = sss_max_tmp;
       sss_idx = i;
   end
   
   sss_max_tmp = sum(xcorr(d5(i,:),subframe)); 
   if(abs(sss_max_tmp) >= sss_max)
       sss_max = sss_max_tmp;
       sss_idx = i;
   end
   
end

sss_idx = sss_idx - 1;

