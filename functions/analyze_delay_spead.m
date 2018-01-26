function [ power_delay ] = analyze_delay_spead(raygain, raydelay, cluster_num, ray_num, Nfft, t_resolution)
%ANALYZE_DELAY_SPEAD Summary of this function goes here
%   Detailed explanation goes here
power_delay = zeros(Nfft,1);

% t_resolution = 4e-9;
hdesign  = fdesign.pulseshaping(t_resolution/1e-9,'Square Root Raised Cosine','Nsym,Beta',7,0.5);
hpulse = design(hdesign);
hpulse_len = length(hpulse.Numerator);
hpulse_ctr = (hpulse_len+1)/2;

for dd=1:Nfft
    for cluster_index = 1:cluster_num
        for ray_index = 1:ray_num
            ofst = abs((raydelay(cluster_index, ray_index)/1e-9) - dd);
            if ofst < hpulse_ctr
               power_delay(dd) = power_delay(dd) + raygain(cluster_index, ray_index)...
                   *hpulse.Numerator(hpulse_ctr + fix(ofst));
            end
        end
    end
end

end

