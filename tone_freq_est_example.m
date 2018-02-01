clear;clc;close all
fc = 1e3;
Ts = 1e-6;
N = 1e3;
SNR_range = 20:5:60;
MCtimes = 1e3;

%%
sigma_fc_range_dB = 0.5:0.5:2; 
sigma_fc_range = [0, 10.^(-sigma_fc_range_dB)];
for sigma_fc_index = 1:length(sigma_fc_range)
    sigma_fc = sigma_fc_range(sigma_fc_index);
    fc_error = zeros(MCtimes, length(SNR_range));
    for MCindex = 1:MCtimes
        for SNRindex = 1:length(SNR_range)
            awgn_pow = 1/(10^(SNR_range(SNRindex)/10));
            n = (1:N)';
            awgn = randn(N,1)*sqrt(awgn_pow);
            for rr=1:20
                dfc = randn*sigma_fc;
                x(:,rr) = exp(1j*Ts*2*pi*(fc+dfc)*(n-1))/20;
            end
            fc_hat = tone_freq_est(sum(x,2)+awgn,Ts);
            fc_error(MCindex, SNRindex) = fc_hat-fc;
        end
    end
    figure(1)
    semilogy(SNR_range,(sqrt(mean(abs(fc_error).^2,1))),'linewidth',2);hold on
end
xlabel('SNR (dB)')
ylabel('RMSE (Hz)')
grid on
title('w/ unknown freq deviation and awgn')
legend('\sigma_f = 0','\sigma_f = 1e-0.5','\sigma_f = 1e-1','\sigma_f = 1e-1.5','\sigma_f = 1e-2')
