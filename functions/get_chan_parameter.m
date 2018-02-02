function [ raygain, raydelay, rayAOA, rayAOD ] = get_chan_parameter( plot_ellipse, print_stat, cluster_num, ray_num)
%GET_CHAN_PARAMETER Summary of this function goes here
%   Detailed explanation goes here

if cluster_num>2
    fprintf('This function only support up to 2 multipath clusters\n');
end

speed_c = 3e8;
D_bs2ue = 40; % Geometric parameter for channel
loc0_ue = [-D_bs2ue/2,0]; % Geometric parameter for channel
loc0_bs = [D_bs2ue/2,0]; % Geometric parameter for channel
cluster_delay = [300e-9,250e-9]; % Mean delay of two multipath clusters


raygain = zeros(cluster_num, ray_num);
raydelay = zeros(cluster_num, ray_num);
rayAOA = zeros(cluster_num, ray_num);
rayAOD = zeros(cluster_num, ray_num);

for cluster_index = 1:cluster_num
    delay = cluster_delay(cluster_index);
    ellipse_len = delay*speed_c;
    a = ellipse_len/2;
    c = D_bs2ue/2;
    b = sqrt(a^2-c^2);
    
%     AOA = rand*pi-pi/2;
    if cluster_index==1
        AOA = pi/2;
    else
        AOA = -pi/2;
    end
    loc_cluster_cnt = [a*cos(AOA),b*sin(AOA)];
    loc_cluster = repmat(loc_cluster_cnt,ray_num,1)+[randn(ray_num,1)*2,randn(ray_num,1)*2];
    loc_cluster_total((cluster_index-1)*ray_num+1:cluster_index*ray_num,:) =  loc_cluster;
    
    for ray_index = 1:ray_num
        raydelay(cluster_index,ray_index) = (norm(loc_cluster(ray_index,:)-loc0_bs)...
                    +norm(loc_cluster(ray_index,:)-loc0_ue))/speed_c;
        temp = loc_cluster(ray_index,:) - loc0_ue;
        rayAOA(cluster_index,ray_index) = angle(temp(1)+1j*temp(2));
        temp = loc_cluster(ray_index,:) - loc0_bs;
        rayAOD(cluster_index,ray_index) = angle(temp(1)+1j*temp(2));
    end
    raygain = exp(1j*rand(cluster_num, ray_num)*2*pi);
%     raygain(2,:) = zeros(1,ray_num);
    
    if print_stat
        fprintf('Cluster %d:\n',cluster_index);
        fprintf('DS Mean:    %4.2f ns\n', mean(raydelay(cluster_index,:)*1e9));
        fprintf('DS Std Dev: %4.2f ns\n', sqrt(var(raydelay(cluster_index,:)*1e9)));
        fprintf('AOAS Mean:    %4.2f degree\n', mean(rayAOA(cluster_index,:)/pi*180));
        fprintf('AOAS Std Dev: %4.2f degree\n', sqrt(var(rayAOA(cluster_index,:)/pi*180)));
        fprintf('AODS Mean:    %4.2f degree\n', mean(rayAOD(cluster_index,:)/pi*180));
        fprintf('AODS Std Dev: %4.2f degree\n', sqrt(var(rayAOD(cluster_index,:)/pi*180)));
    end
    
    if plot_ellipse % plot ellipse
        t = -pi:0.01:pi;
        x = a*cos(t);
        y = b*sin(t);
        figure(99)
        plot(x,y,'g--');hold on
    end
end

if plot_ellipse % plot of scattering ellipse
    figure(99)
    plot(loc0_ue(1),loc0_ue(2),'x');hold on
    plot(loc0_bs(1),loc0_bs(2),'x');hold on
    plot(loc_cluster_total(:,1),loc_cluster_total(:,2),'o');hold on
    title('Geometric Stochastic Multipath')
    axis([-80,80,-80,80])
    xlabel('X Coordn (m)')
    ylabel('Y Coordn (m)')
    grid on
end


end

