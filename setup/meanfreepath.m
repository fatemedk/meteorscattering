% plot mean free path versus altitude

clear all; close all;

loadconstants;

alt = 70:1:140;

% get neutral density

nd = MSISatmosphere1(alt);
T = nd.temp;

% need to know viscosity: use Sutherland's formula

C = 120;
T0 = (5/9)*524.07;
mu0 = 1.827e-5;
a = T0 + C;
b = T + C;
mu = mu0*(a./b).*(nd.temp/T0).^(1.5);

figure; plot(nd.temp,alt,mu*1e7,alt);
legend('Temperature (K)','Viscosity * 1e7');

p = nd.total*1e6 * kB .* nd.temp;
p0 = 1.0103e5;

% use relative-to-ground values, and assume mean free path at sea level is 65 nm

L0 = 65e-9;

L = L0 * (mu/mu0) .* (p0./p) .* sqrt(T/T0);

figure; semilogx(L,alt);
grid on;
hold on;

% plot wavelengths of different radars

radarfreq = [50e6 160e6 422e6 450e6];
radarlam = 3e8./radarfreq;
colors = colormap(jet(length(radarlam))) * 0.8;

for m = 1:length(radarlam),
    plot([radarlam(m) radarlam(m)],[alt(1) alt(end)],'color',colors(m,:));
end

legend('Mean Free Path (m)','JRO','Altair VHF','Altair UHF','PFISR');