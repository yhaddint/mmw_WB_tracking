function [ loc0_ue, loc0_bs, loc_cluster_total ] = get_init_locations( plot_ellipse, cluster_num, ray_num, scatter_radius )
%GET_LOCOCATIONS Summary of this function goes here
%   Detailed explanation goes here
    if nargin==3
        scatter_radius = 1.5;
    end
    
    speed_c = 3e8;
    D_bs2ue = 50; % Geometric parameter for channel
    loc0_ue = [-D_bs2ue/2,0]; % Geometric parameter for channel
    loc0_bs = [D_bs2ue/2,0]; % Geometric parameter for channel
    cluster_delay = [200e-9,250e-9]; % Mean delay of two multipath clusters

    
    for cluster_index = 1:cluster_num
        delay = cluster_delay(cluster_index);
        ellipse_len = delay*speed_c;
        a = ellipse_len/2;
        c = D_bs2ue/2;
        b = sqrt(a^2-c^2);

    %     AOA = rand*pi-pi/2;
        if cluster_index == 1
            AOA = 0.6*pi;
        else
            AOA = -pi/2;
        end
        loc_cluster_cnt = [a*cos(AOA),b*sin(AOA)];
        loc_cluster = repmat(loc_cluster_cnt,ray_num,1)+[randn(ray_num,1)*scatter_radius,randn(ray_num,1)*scatter_radius];
        loc_cluster_total((cluster_index-1)*ray_num+1:cluster_index*ray_num,:) =  loc_cluster;


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

