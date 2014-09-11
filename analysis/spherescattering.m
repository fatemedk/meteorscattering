% sphere scattering: compare to theory

clear all; close all;

topdir = '/shared/users/ram80/meteor/runs/sphere/';
doplot = 0;

loadconstants;

sizes = [0.05 0.075 0.1 0.12 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.8 1 1.2 1.5 1.8 2 2.5 3 4];

h = figure(1);
ax = axes;

for m = 1:length(sizes),
    
    datadir = sprintf('%splatesize_%g/',topdir,sizes(m));
    load([datadir 'inputs.mat']);
    
    if exist([datadir 'rcs.mat'],'file'),
        load([datadir 'rcs.mat']);
    else
        sig = meteorNTF(datadir,doplot);
    end
    
    Clam = 2*pi*in.platesize./(vp./in.fntf);
    loglog(Clam,sig.sig./(pi*in.platesize^2),'rx-');
    hold(ax,'on');
    sig.expected = sphereRCSvsFreq(in.fntf,in.platesize);
    loglog(Clam,sig.expected./(pi*in.platesize^2),'kx-');
    
    drawnow;
end

%%

set(ax,'xlim',[1e-1 1e1],'ylim',[1e-3 1e1]);

xlabel(ax,'Sphere Circumference / lambda');
ylabel(ax,'RCS / (pi * r^2)');
title(ax,'Radar Cross Section for conducting sphere');
legend(ax,'FDTD Method','Mie theory');
