% Quick look at main-lobe and side-lobes of spatial beams

clear;clc;
Nr = 32;
phi0 = 30/180*pi;
phi_range_deg = (-90:0.1:90);
phi_range_rad = phi_range_deg/180*pi;
for pp=1:length(phi_range_rad)
    phi = phi_range_rad(pp);
    bigPhi = sin(phi0)-sin(phi);
    bigPhi_aprx = phi0-phi;
    BP(pp) = (1-exp(1j*pi*Nr*bigPhi))/(1-exp(1j*pi*bigPhi));
    BP2(pp) = (1-exp(1j*pi*Nr*bigPhi_aprx))/(1-exp(1j*pi*bigPhi_aprx));

end
figure;
plot(phi_range_deg,20*log10(abs(BP)));hold on
plot(phi_range_deg,20*log10(abs(BP2)));hold on
grid on
ylim([-10,40])
xlim([-90,90])
