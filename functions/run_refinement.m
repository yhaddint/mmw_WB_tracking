function [ phi_last ] = run_refinement( chan_noisy_ob, last_rayAOA, last_rayAOD, F, W, cluster_num, Nt, Nr )
%RUN_REFINEMENT Summary of this function goes here
%   Detailed explanation goes here
    ite_num = 25;
    
    phi_hat(1) = last_rayAOA;
    
    for ii = 1:ite_num
        %---------------------------------------------
        % Alpha estimation using previous Phi
        %---------------------------------------------
        arx_hat = exp(1j*(0:Nr-1)'*pi*sin(phi_hat(ii)))/sqrt(Nr);
        atx_hat = exp(1j*(0:Nt-1)'*pi*sin(last_rayAOD))/sqrt(Nt);

        H_cal = diag(W' * arx_hat * atx_hat' * F);
        alpha_hat(ii) = pinv(H_cal) * chan_noisy_ob;
        error(ii) = norm(alpha_hat(ii)*H_cal - chan_noisy_ob);
        %---------------------------------------------
        % Phi estimation using previous Alpha
        %---------------------------------------------      
        dPhi = exp(1j*pi*(0:Nr-1).'*sin(phi_hat(ii)))/sqrt(Nr).*(1j*pi*(0:Nr-1).'*cos(phi_hat(ii)));
        
        D_cal = alpha_hat(ii) * diag(W' * dPhi * atx_hat' * F);
        yr = chan_noisy_ob - H_cal * alpha_hat(ii);
        Q_cal = [real(D_cal);imag(D_cal)];
        tr = [real(yr);imag(yr)];

        dphi = pinv(Q_cal) * tr;
        phi_hat(ii+1) = phi_hat(ii) + dphi;
        
    end
    phi_last = phi_hat(end);
end

