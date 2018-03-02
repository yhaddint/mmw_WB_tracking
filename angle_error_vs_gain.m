clear;clc;
ray_num = 4;
sigma_AOA_spread = 10/180*pi;
cluster_index = 1;
sweep_num = 100;
sweep_range = linspace(-10,10,sweep_num)/180*pi;
MCtimes = 500;
Nt = 16;

for MCindex = 1:MCtimes
    ray_AOA_azim_mean = 0;
    angle_spread = laprnd(1, ray_num, 0, sigma_AOA_spread);
    ray_AOA_azim(cluster_index,:) = ray_AOA_azim_mean + angle_spread;

    chan_MISO = zeros(Nt,1);
    for rr = 1:ray_num
%         g = (randn + 1j * randn)/sqrt(ray_num);
        g = 1/ray_num;

        chan_MISO = chan_MISO + 1/ray_num*exp(1j*pi*(0:Nt-1).'*sin(ray_AOA_azim(rr)));
    end
    chan_MISO_norm = chan_MISO./norm(chan_MISO)*sqrt(Nt);
    
    BF_opt = chan_MISO_norm./abs(chan_MISO_norm);
    result_opt(MCindex) = abs(chan_MISO_norm'*BF_opt)^2/Nt;
    
    for kk=1:sweep_num
        steer_vec = exp(1j*pi*(0:Nt-1).'*sin(sweep_range(kk)));
        result(kk,MCindex) = abs(chan_MISO_norm'*steer_vec)^2/Nt;
    end
end

figure
plot(sweep_range/pi*180,10*log10(mean(result,2)));hold on
plot(sweep_range/pi*180,ones(sweep_num,1)*10*log10(mean(result_opt)));hold on
grid on
xlabel('Angle Error for Centroid (deg)')
ylabel('Gain (dB)')
legend('Steering at Centroid Est.','Steer based on True CSI')