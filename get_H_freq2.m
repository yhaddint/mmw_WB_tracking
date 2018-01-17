function [ H_freq ] = get_H_freq2( raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr )
%GET_H_FREQ Summary of this function goes here
%   Detailed explanation goes here

    Nfft = 512;
    H_freq = zeros(Nr,Nt,Nfft);
    for cluster_index = 1:cluster_num
        for ray_index = 1:ray_num
            delaynow = raydelay(cluster_index, ray_index);
            phi = rayAOA(cluster_index, ray_index);
            arx = exp(1j*(0:Nr-1)'*pi*sin(phi))/sqrt(Nr);
            theta = rayAOD(cluster_index, ray_index);
            atx = exp(1j*(0:Nt-1)'*pi*sin(theta))/sqrt(Nt);
            for kk=1:Nfft
            H_freq(:,:,kk) = H_freq(:,:,kk) ...
                + exp(-1j*2*pi*kk/Nfft*delaynow/(1e-9))*raygain(cluster_index, ray_index)*arx*atx';
            end

        end
    end

end