%% setup file: change the parameters as seen fit, then run this script 
% to create the inputs.dat file that will be read by emp2d. 

function [in,jobid] = meteorrun(in)

loadconstants;

in.metsize = [in.metsize in.metsize*in.metshape];

in.lambda = 3e8/in.f0;

% change gridfactor so that dx is at least metsize/50

%tempdx = in.lambda/in.gridfactor;
%if (tempdx > in.metsize(1)/50),
%    newdx = in.metsize(1)/50;
%    in.gridfactor = in.lambda/newdx;
%end

% finalized grid cell size

in.dx = in.lambda/in.gridfactor;
in.dy = in.dx; 
in.dz = in.dx;

% grid cell locations
in.x = -(in.Nx-1)/2*in.dx + in.dx*(0:in.Nx-1);
in.y = -(in.Ny-1)/2*in.dy + in.dy*(0:in.Ny-1);
in.z = -(in.Nz-1)/2*in.dz + in.dz*(0:in.Nz-1);

% calculate time steps from other parameters

in.stabfactor = 1 / sqrt(1 + (in.wpmax*pi/in.gridfactor/2)^2);
in.dt = in.dx/vp/sqrt(3) * in.stabfactor * 0.8;
in.tsteps = round(in.Nx*in.dx/vp/in.dt*in.tfactor);


%% create input file

if ~exist(in.rundir,'dir'),
    mkdir(in.rundir);
end

fid = fopen([in.rundir '/inputs.dat'],'w');
fwrite(fid,in.dopml,'int');
fwrite(fid,in.savefields,'int');
fwrite(fid,in.f0,'double');
fwrite(fid,in.fsig,'double');
fwrite(fid,in.sourcetype,'int');
fwrite(fid,in.sourcepol,'int');
fwrite(fid,in.gridfactor,'double');
fwrite(fid,in.dt,'double');
fwrite(fid,in.tsteps,'int');
fwrite(fid,[in.Nx in.Ny in.Nz],'int');
fwrite(fid,in.metsize,'double');
fwrite(fid,in.wpmax,'double');
fwrite(fid,in.numax,'double');
fwrite(fid,in.metloc,'double');
fwrite(fid,in.metangle,'double');
fwrite(fid,in.mettype,'int');
fwrite(fid,in.numfiles,'int');
fwrite(fid,in.nprobes,'int');
fwrite(fid,in.probex,'int');
fwrite(fid,in.probey,'int');
fwrite(fid,in.probez,'int');
fwrite(fid,in.nntf,'int');
fwrite(fid,in.fntf,'double');
fwrite(fid,in.doplate,'int');
fwrite(fid,in.platesize,'double');
fclose(fid);

% magnetic field file

fid = fopen([in.rundir '/B0.dat'],'w');
fwrite(fid,in.Bx,'double');
fwrite(fid,in.By,'double');
fwrite(fid,in.Bz,'double');
fclose(fid);


%% create meteor plasma

wp = zeros(in.Nx,in.Ny);
nu = zeros(in.Nx,in.Ny);

msig = in.metsize;

for i = 1:in.Nx,
    for j = 1:in.Ny,
        if in.x(i) < in.metloc(1),
        wp(i,j) = in.wpmax*(2*pi*in.f0) * exp(-(in.x(i) - in.metloc(1))^2/msig(1)^2) * exp(-(in.y(j) - in.metloc(2))^2/msig(1)^2);
        else
        wp(i,j) = in.wpmax*(2*pi*in.f0) * exp(-(in.x(i) - in.metloc(1))^2/msig(2)^2) * exp(-(in.y(j) - in.metloc(2))^2/msig(1)^2);    
        end
    end
end

% rotate by theta!
in.wp = imrotate(wp,in.metangle,'crop');

% get collision frequencies
in = meteorCollisions(in);


% add electron-neutral collision frequency to coulomb collisions above
in.nu = in.coll.en + in.coll.ee + in.coll.ei;

fid = fopen([in.rundir '/wp.dat'],'w');
fwrite(fid,in.wp,'double');
fclose(fid);

fid = fopen([in.rundir '/nu.dat'],'w');
fwrite(fid,in.nu,'double');
fclose(fid);


% plot them

h1 = figure(1);
set(h1,'position',[200 400 1800 600]);
ax1 = subplot(131);
ax2 = subplot(132);
ax3 = subplot(133);
imagesc(in.x,in.y,in.wp,'parent',ax1); axis xy; colorbar('peer',ax1);
imagesc(in.x,in.y,in.nu,'parent',ax2); axis xy; colorbar('peer',ax2);

hold(ax1,'on');

contour(ax1,in.x,in.y,in.wp'/(2*pi*in.f0),[1, 1],'Color',[1 1 1]);
contour(ax1,in.x,in.y,in.wp'/(2*pi*in.f0),in.wpmax*exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);

% plot source

[source,ax3] = plotMeteorSource(in,ax3);


%%  measure wp/w = 1 radius in y-dimension

wpslice = max(in.wp,[],1);

i1 = find(wpslice > (2*pi*in.f0),1,'first');
i2 = find(wpslice > (2*pi*in.f0),1,'last');

in.wpdiameter = (i2-i1)*in.dy;


%% save inputs to mat file

save([in.rundir '/inputs.mat'],'in');


%% run the simulation

if (in.submitjob),
    
    % create pbs file to run simulation
    
    pbsfile = writepbsfile(in.rundir,in.runname,in.exefile);
    
    % run command
    
    system(['cp ' in.exedir in.exefile ' ' in.rundir]);
    
    % for simplicity, cd into run directory, run it, then return to pwd
    
    thisdir = pwd;
    cd(in.rundir);
    
    submitstr = ['qsub -q ' in.cluster ' -d ' in.rundir ' -l nodes=1:ppn=' in.numnodes ' -l walltime=72:00:00 ' ...
        pbsfile];
    
    [~,jobname] = system(submitstr);
    jobid = strtrim(jobname);
    
    fprintf('Job %s submitted!\n',jobid);
    
    cd(thisdir);
    
else
    
    jobid = '';
    fprintf('Job not submitted\n');
    
end
