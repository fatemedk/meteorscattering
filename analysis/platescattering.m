% plot of radar cross sections for different plate scatterers; to show off
% how my NTF is working.

clear all; close all;

plate = [1 2 4 6 8 10 12 14];

load('/shared/users/ram80/meteor/source/inputs.mat');

h1 = figure(1);
ax = axes;
colors = colormap(jet(length(plate))) * 0.8;

% note dBsm is 10*log10!

p = zeros(2*length(plate),1);

for m = 1:length(plate),
    
    load(sprintf('sig_plate_%02.0f',plate(m)));
    
    p(2*m-1) = plot(ax,in.fntf/1e6,10*log10(sig.sig),'o-','color',colors(m,:));
    hold on;
    p(2*m) = plot(ax,in.fntf/1e6,10*log10(sig.expected),'x-','color',colors(m,:));
    
end

xlabel(ax,'Frequency (MHz)','Fontsize',16);
ylabel(ax,'Radar cross section (dBsm)','Fontsize',16);
title(ax,'Scattering Cross Section for different flat plates','Fontsize',16);
    
legend(ax,p([1 2 3:2:end]),'1x1 m plate','Theoretical','2x2 m','4x4 m','6x6 m','8x8 m','10x10 m','12x12 m','14x14 m');

set(ax,'Fontsize',16);