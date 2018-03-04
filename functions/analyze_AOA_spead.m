function [ power_AOA ] = analyze_AOA_spead( H_freq, Nt, Nr, Nfft)
%ANALYZE_AOA_SPEAD Summary of this function goes here
%   Detailed explanation goes here
    if nargin==3
        power_AOA = zeros(32,1);
        phi_range = linspace(-pi/2,pi/2,256);
        for pp=1:256
            phi = phi_range(pp);
            arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
            power_subcarrier = norm(arx'*H_freq)^2;
            power_AOA(pp) = power_subcarrier;
        end
    else
        phi_range = linspace(-pi/2,pi/2,64);
        for pp=1:64
            phi = phi_range(pp);
            arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
            for kk=1:Nfft
                power_subcarrier(kk) = norm(arx'*squeeze(H_freq(:,:,kk)))^2;
            end
            power_AOA(pp) = sum(power_subcarrier);
        end
    end
end
