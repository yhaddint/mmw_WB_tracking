function [ tx_max_index, tx_angle, rx_max_index, rx_angle ] = run_angle_sweep_2D( Nt, Nr, Nfft, chan_WB, angle_tx_range, angle_rx_range)
%RUN_ANGLE_SWEEP_2D Summary of this function goes here
%   Detailed explanation goes here

    rx_RSS_bin = zeros(Nfft,1);
    rx_RSS = zeros(length(angle_rx_range),length(angle_tx_range));

    for beam_tx_index = 1:length(angle_tx_range)
        theta = angle_tx_range(beam_tx_index);
        atx = exp(1j*(0:Nt-1)'*pi*sin(theta))/sqrt(Nt);
        for beam_rx_index = 1:length(angle_rx_range)
            phi = angle_rx_range(beam_rx_index);
            arx = exp(1j*(0:Nr-1)'*pi*sin(phi))/sqrt(Nr);
            for kk=1:Nfft
                rx_RSS_bin(kk) = arx'*(squeeze(chan_WB(:,:,kk))) * atx;
            end
            rx_RSS(beam_rx_index,beam_tx_index) = 10*log10(sum(abs(rx_RSS_bin).^2));
        end
    end
    [~,tx_max_index] = max(max(rx_RSS));
    [~,rx_max_index] = max(max(transpose(rx_RSS)));
    
    tx_angle = angle_tx_range(tx_max_index)/pi*180;
    rx_angle = angle_rx_range(rx_max_index)/pi*180;
end