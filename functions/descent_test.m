function [ii] = descent_test( chan_noisy_ob, last_rayAOA, last_rayAOD, F, W, cluster_num, Nt, Nr )
%DESCENT_TEST Summary of this function goes here
%   Detailed explanation goes here
    phi_range = linspace(-5,5,100)/180*pi;
    
    for ii = 1:length(phi_range)
        %---------------------------------------------
        % Alpha estimation using previous Phi
        %---------------------------------------------
        arx_hat = exp(1j*(0:Nr-1)'*pi*sin(phi_range(ii)))/sqrt(Nr);
        atx_hat = exp(1j*(0:Nt-1)'*pi*sin(last_rayAOD))/sqrt(Nt);

        H_cal = diag(W' * arx_hat * atx_hat' * F);
        alpha_hat(ii) = pinv(H_cal) * chan_noisy_ob;
        error(ii) = norm(alpha_hat(ii)*H_cal - chan_noisy_ob);
    end
    figure
    plot(phi_range/pi*180,sqrt(error+1e-10)/pi*180);
    grid on

end

