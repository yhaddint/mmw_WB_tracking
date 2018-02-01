function [ fc_hat ] = tone_freq_est( sig, Ts )
%TONE_FREQ_EST Summary of this function goes here
%   Detailed explanation goes here

    L = length(sig);
    phase_est = mean(phase(sig(2:L) .* conj(sig(1:L-1))));
    fc_hat = phase_est/(Ts*2*pi);
end

