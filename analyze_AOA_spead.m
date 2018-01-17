function [ power_AOA ] = analyze_AOA_spead( H_freq, Nt, Nr, Nfft)
%ANALYZE_AOA_SPEAD Summary of this function goes here
%   Detailed explanation goes here
phi_range = linspace(-pi/2,pi/2,64);
for pp=1:64
    phi = phi_range(pp);
    arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
    for kk=1:Nfft
        power_subcarrier(kk) = norm(arx'*squeeze(H_freq(:,:,kk)))^2;
    end
    power_AOA(pp) = sum(power_subcarrier);
end

