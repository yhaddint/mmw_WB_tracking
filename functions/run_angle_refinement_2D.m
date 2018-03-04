function [ current_rx_opt_index ] = run_angle_refinement_2D( Nr, prev_rx_opt_index, atx_opt, Nfft, chan_WB, angle_rx_range,noise_power)
%RUN_ANGLE_REFINEMENT_2D Summary of this function goes here
%   Detailed explanation goes here
    if nargin==6
        noise_power = 0;
    end
    if  Nfft==1
        refine_index_range = [0, 1, -1, 2, -2];
        refine_num = length(refine_index_range);
        for refine_index = 1:refine_num
            rx_refine_index = prev_rx_opt_index + refine_index_range(refine_index);
%             if rx_refine_index>16
%                 rx_refine_index=16;
%             end
            phi_refine = angle_rx_range(rx_refine_index);
            arx_refine = exp(1j*(0:Nr-1)'*pi*sin(phi_refine));
            noise = randn + 1j*randn;
            awgn = noise * sqrt(noise_power/2);
            rx_RSS_refine_bin = arx_refine'* chan_WB * atx_opt + awgn;
            rx_refine_RSS(refine_index) = 10*log10(abs(rx_RSS_refine_bin).^2);
        end
        [~, opt_index] = max(rx_refine_RSS);
        current_rx_opt_index = prev_rx_opt_index + refine_index_range(opt_index);
    else
        refine_index_range = [0, 1, -1, 2, -2];
        refine_num = length(refine_index_range);
        for refine_index = 1:refine_num
            rx_refine_index = prev_rx_opt_index + refine_index_range(refine_index);
            phi_refine = angle_rx_range(rx_refine_index);
            arx_refine = exp(1j*(0:Nr-1)'*pi*sin(phi_refine));
            for kk=1:Nfft
                rx_RSS_refine_bin(kk) = arx_refine'*(squeeze(chan_WB(:,:,kk)))*atx_opt;
            end
            rx_refine_RSS(refine_index) = 10*log10(sum(abs(rx_RSS_refine_bin).^2));
        end
        [~, opt_index] = max(rx_refine_RSS);
        current_rx_opt_index = prev_rx_opt_index + refine_index_range(opt_index);
    end
end

