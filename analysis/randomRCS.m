% create scatter plot of RCS versus size, wpmax for random meteor
% simulations

clear all; close all;

topdir = '/shared/users/ram80/meteor/runs/randomset/';

d = dir([topdir 'random*']);

for m = 1:length(d),
    
    fprintf('Doing directory %s...\n',d(m).name);
    
    datadir = [topdir d(m).name '/'];
    load([datadir 'inputs.mat']);
    if exist([datadir 'rcs.mat'],'file'),
        load([datadir 'rcs.mat']);
    else
        sig = meteorNTF(datadir,0);
    end
    
    msize(m) = in.metsize(1);
    wpmax(m) = in.wpmax;
    n = length(sig.sig);
    rcs(m) = sig.sig((n+1)/2);
    
    fprintf('Size = %.3f, wpmax = %.3f, sig = %.3f\n',...
        msize(m),wpmax(m),rcs(m));
    
end

%% 2D scatter plot, with markers by color

h1 = figure(1);

ax1 = subplot(221);
scatter(ax1,msize,wpmax,100,20*log10(rcs),'filled');
xlabel(ax1,'Meteor radius (m)');
ylabel(ax1,'Peak wp/w0');
cax1 = colorbar('peer',ax1);
ylabel(cax1,'RCS (dBsm)');

sizefactor = wpmax.*msize.^2;

ax2 = subplot(222);
scatter(ax2,msize,20*log10(rcs));
xlabel(ax2,'Meteor size (m)');
ylabel(ax2,'RCS (dBsm)');

ax3 = subplot(223);
scatter(ax3,wpmax,20*log10(rcs));
xlabel(ax3,'Max wp/w0');
ylabel(ax3,'RCS 9dBsm)');

ax4 = subplot(224);
scatter(ax4,log10(sizefactor),20*log10(rcs));
xlabel(ax4,'Log10 Size factor');
ylabel(ax4,'RCS 9dBsm)');
title(ax4,'Size factor = wpmax * size^2');

