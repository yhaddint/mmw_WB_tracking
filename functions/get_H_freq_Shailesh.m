function [ H_freq, Delta_freq] = get_H_freq_Shailesh(raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr, Nfft, Ts )                                   
                                            
%GET_H_FREQ Summary of this function goes here
%   Detailed explanation goes here
    H_time = zeros(Nr,Nt,Nfft);
    Delta_D = zeros(1,Nfft);
    for dd=1:Nfft
        for cluster_index = 1:cluster_num
            for ray_index = 1:ray_num
                if abs(raydelay(cluster_index, ray_index)-(dd*Ts))<Ts/2
                    phi = rayAOA(cluster_index, ray_index);
                    arx = exp(1j*(0:Nr-1)'*pi*sin(phi));
                    theta = rayAOD(cluster_index, ray_index);
                    atx = exp(1j*(0:Nt-1)'*pi*sin(theta));
                    H_time(:,:,dd) = H_time(:,:,dd) + arx*atx';
                    Delta_D(dd) = 1; % 1 because it is the multiplying factor to arx*atx';
                end
            end
        end
    end

    
    H_freq = zeros(Nr,Nt,Nfft);
    Delta_freq = fft(Delta_D);
    for n1=1:Nr
        for n2=1:Nt
            H_SISO_time = squeeze(H_time(n1,n2,:));
            H_freq(n1,n2,:) = fft(H_SISO_time);
        end
    end

%  k = 10;
%  H_freq(:,:,k)
%  arx*Delta_freq(k)*atx';
%  
% 
% a=1;

