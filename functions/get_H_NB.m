function [ H_NB ] = get_H_NB( raygain, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr )
%GET_H_FREQ Summary of this function goes here
%   This is the function used for most script
%   Channel model is H = sum
%   alpha*exp(-j*2*pi*n*tau_0/(Nfft*Ts))*ar(phi)*at(theta)
%   H_freq = get_H_freq2( raygain, raydelay, rayAOA, rayAOD, cluster_num, ray_num, Nt, Nr )
%   IP: raygain (cluster_num by ray_num) matrix with complex gain for each path
%   IP: raydelay (cluster_num by ray_num) matrix with prop delay for each path
%   IP: rayAOA (cluster_num by ray_num) matrix with AOA for each path
%   IP: rayAOD (cluster_num by ray_num) matrix with AOD for each path
%   IP: cluster_num is scaler for number of multipath cluster
%   IP: ray_num is scaler for number of rays in each cluster
%   IP: Nt is number of antenna (ULA) in transmitter
%   IP: Nr is number of antenna (ULA) in receiver

    H_NB = zeros(Nr,Nt);
    for cluster_index = 1:cluster_num
        for ray_index = 1:ray_num
            phi = rayAOA(cluster_index, ray_index);
            arx = exp(1j*(0:Nr-1)'*pi*sin(phi))/sqrt(Nr);
            theta = rayAOD(cluster_index, ray_index);
            atx = exp(1j*(0:Nt-1)'*pi*sin(theta))/sqrt(Nt);           
            H_NB = H_NB + raygain(cluster_index, ray_index)*arx*atx';

        end
    end

end