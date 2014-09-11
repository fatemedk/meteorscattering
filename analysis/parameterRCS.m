% create scatter plot of RCS versus size, wpmax parameter sweeps

clear all; close all;

topdir = '/shared/users/ram80/meteor/runs/wpmax/';

d = dir([topdir 'wpmax*']);

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

% scatter plot, with markers by color

h1 = figure(1);

ax1 = subplot(121);
plot(ax1,wpmax,20*log10(rcs),'x-');
xlabel(ax1,'Peak wp/w0');
ylabel(ax1,'RCS (dBsm)');

%% do it again for meteor size

topdir = '/shared/users/ram80/meteor/runs/size/';

d = dir([topdir 'metsize*']);

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

ax2 = subplot(122);
plot(ax2,msize,20*log10(rcs),'x-');
xlabel(ax2,'Meteor size (m)');
ylabel(ax2,'RCS (dBsm)');


