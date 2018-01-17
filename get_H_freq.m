function [ H_freq ] = get_H_freq( raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr, fc )
%GET_H_FREQ Summary of this function goes here
%   Detailed explanation goes here
    H_time = zeros(Nr,Nt,512);
    for dd=1:512
        for cluster_index = 1:cluster_num
            for ray_index = 1:ray_num
                delaynow = raydelay(cluster_index, ray_index);
                if abs(delaynow-(dd*1e-9))<0.5e-9
                    phi = rayAOA(cluster_index, ray_index);
                    arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
                    theta = rayAOD(cluster_index, ray_index);
                    atx = exp(1j*(0:Nt-1)'*pi*sin(theta));
%                     H_time(:,:,dd) = H_time(:,:,dd) + exp(-1j*fc*delaynow*2*pi)*arx*atx';
                    H_time(:,:,dd) = H_time(:,:,dd) ...
                        + raygain(cluster_index, ray_index)*arx*atx';

                end
            end
        end
    end

    Nfft = 512;
    H_freq = zeros(Nr,Nt,Nfft);
    for n1=1:Nr
        for n2=1:Nt
            H_SISO_time = squeeze(H_time(n1,n2,:));
            H_freq(n1,n2,:) = ifft(H_SISO_time);
        end
    end

end

