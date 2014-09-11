% quickly read results from my 3D code

clear all; close all;

runname = 'metsize_0.5';
datadir = ['/shared/users/ram80/meteor/runs/size/' runname '/'];

saveit = 1;

% choose either 'total' or 'scattered' for second row of plot fields. If
% total, plots Z component (first row is Y); if scattered, plots Y
% scattered.

plotfield2 = 'scattered';

loadconstants;

load([datadir 'inputs.mat']);

dx = in.dx;
dy = in.dy;
dz = in.dz;
x = linspace(-dx*in.Nx/2,dx*in.Nx/2,in.Nx);
y = linspace(-dy*in.Ny/2,dy*in.Ny/2,in.Ny);
z = linspace(-dz*in.Nz/2,dz*in.Nz/2,in.Nz);
w0 = 2*pi*in.f0;
% meteor size in grid cells
metsize = in.metsize / in.dx;
doplate = in.doplate;

in.tsteps = round(in.Nx*in.dx/vp/in.dt*3.0);
writedt = (in.tsteps/in.numfiles)*in.dt;

% tf/sf boundaries
bx1 = 21; bx2 = in.Nx-bx1+1;
by1 = 21; by2 = in.Ny-by1+1;
bz1 = 21; bz2 = in.Nz-bz1+1;

% I want to plot meteor lines; for that I need to draw meteor
wpxy = zeros(in.Nx,in.Ny);
wpxz = wpxy;
wpyz = wpxy;
for i = 1:in.Nx,
    for j = 1:in.Ny,
        if i < in.Nx/2,
        wpxz(i,j) = (in.wpmax *w0) * exp(-(i - (in.Nx+1)/2)^2/metsize(1)^2) * exp(-(j - (in.Ny+1)/2)^2/metsize(1)^2);
        else
        wpxz(i,j) = (in.wpmax *w0) * exp(-(i - (in.Nx+1)/2)^2/metsize(2)^2) * exp(-(j - (in.Ny+1)/2)^2/metsize(1)^2);    
        end
        wpyz(i,j) = (in.wpmax *w0) * exp(-(i - (in.Nx+1)/2)^2/metsize(1)^2) * exp(-(j - (in.Ny+1)/2)^2/metsize(1)^2);
    end
end

fid = fopen([datadir 'wp.dat'],'r');
wpxy = fread(fid,[in.Nx in.Ny],'double');
fclose(fid);


% figure setup
h1 = figure(1); set(h1,'position',[30 30 1200 900]);
set(h1,'PaperUnits','points','PaperPosition',get(h1,'position'));
set(h1,'color',[1 1 1]);

set(h1,'renderer','zbuffer');

ax1 = subplot(231);
ax2 = subplot(232);
ax3 = subplot(233);
ax4 = subplot(234);
ax5 = subplot(235);
ax6 = subplot(236);

c = colormap('jet');
c3 = [0 0 0; 0 0 0.1; 0 0 0.2; 0 0 0.3; 0 0 0.4; c];
c2 = [1 1 1; 0.8 0.8 0.9; 0.6 0.6 0.8; 0.4 0.4 0.7; 0.2 0.2 0.6; c];

% open files for reading

Efile = fopen([datadir 'output_E.dat'],'r');
Jfile = fopen([datadir 'output_J.dat'],'r');
Hfile = fopen([datadir 'output_H.dat'],'r');

EYslicez = zeros(in.Nx,in.Ny);
EYslicey = zeros(in.Nx,in.Nz);
EYslicex = zeros(in.Ny,in.Nz);
JYslicez = zeros(in.Nx,in.Ny);
JYslicey = zeros(in.Nx,in.Nz);
JYslicex = zeros(in.Ny,in.Nz);
HYslicez = zeros(in.Nx,in.Ny);
HYslicey = zeros(in.Nx,in.Nz);
HYslicex = zeros(in.Ny,in.Nz);

EZslicez = zeros(in.Nx,in.Ny);
EZslicey = zeros(in.Nx,in.Nz);
EZslicex = zeros(in.Ny,in.Nz);
JZslicez = zeros(in.Nx,in.Ny);
JZslicey = zeros(in.Nx,in.Nz);
JZslicex = zeros(in.Ny,in.Nz);
HZslicez = zeros(in.Nx,in.Ny);
HZslicey = zeros(in.Nx,in.Nz);
HZslicex = zeros(in.Ny,in.Nz);

t = 0:(in.tsteps-1);

% initialize plots

im1 = imagesc(x,y,log10(abs(EYslicez)+1e-8),'parent',ax1); %shading(ax1,'flat');
im2 = imagesc(x,z,log10(abs(EYslicey)+1e-8),'parent',ax2); %shading(ax2,'flat');
im3 = imagesc(y,z,log10(abs(EYslicex)+1e-8),'parent',ax3); %shading(ax3,'flat');
im4 = imagesc(x,y,log10(abs(HYslicez)+1e-8),'parent',ax4); %shading(ax4,'flat');
im5 = imagesc(x,z,log10(abs(HYslicey)+1e-8),'parent',ax5); %shading(ax5,'flat');
im6 = imagesc(y,z,log10(abs(HYslicex)+1e-8),'parent',ax6); %shading(ax6,'flat');

set(ax1,'Fontsize',16);
set(ax2,'Fontsize',16);
set(ax3,'Fontsize',16);
set(ax4,'Fontsize',16);
set(ax5,'Fontsize',16);
set(ax6,'Fontsize',16);

ti1 = title(ax1,'','Fontsize',16);
ti2 = title(ax2,'','Fontsize',16);
ti3 = title(ax3,'','Fontsize',16);
ti4 = title(ax4,'','Fontsize',16);
ti5 = title(ax5,'','Fontsize',16);
ti6 = title(ax6,'','Fontsize',16);

%colormap(ax1,c2);

caxis(ax1,[-5 0]);
caxis(ax2,[-5 0]);
caxis(ax3,[-5 0]);
caxis(ax4,[-5 0]);
caxis(ax5,[-5 0]);
caxis(ax6,[-5 0]);

xlabel(ax1,'x (meters)','Fontsize',16); ylabel(ax1,'y (meters)','Fontsize',16);
xlabel(ax2,'x (meters)','Fontsize',16); ylabel(ax2,'z (meters)','Fontsize',16);
xlabel(ax3,'y (meters)','Fontsize',16); ylabel(ax3,'z (meters)','Fontsize',16);
xlabel(ax4,'x (meters)','Fontsize',16); ylabel(ax4,'y (meters)','Fontsize',16);
xlabel(ax5,'x (meters)','Fontsize',16); ylabel(ax5,'z (meters)','Fontsize',16);
xlabel(ax6,'y (meters)','Fontsize',16); ylabel(ax6,'z (meters)','Fontsize',16);

hold(ax1,'on');
hold(ax2,'on');
hold(ax3,'on');
hold(ax4,'on');
hold(ax5,'on');
hold(ax6,'on');

if ~doplate,
    
contour(ax1,x,y,wpxy'/w0,[1, 1],'Color',[1 1 1]);
contour(ax2,x,z,wpxz'/w0,[1, 1],'Color',[1 1 1]);
contour(ax3,y,z,wpyz'/w0,[1, 1],'Color',[1 1 1]);
contour(ax4,x,y,wpxy'/w0,[1, 1],'Color',[1 1 1]);
contour(ax5,x,z,wpxz'/w0,[1, 1],'Color',[1 1 1]);
contour(ax6,y,z,wpyz'/w0,[1, 1],'Color',[1 1 1]);

contour(ax1,x,y,wpxy'/in.wpmax/w0,exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);
contour(ax2,x,z,wpxz'/in.wpmax/w0,exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);
contour(ax3,y,z,wpyz'/in.wpmax/w0,exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);
contour(ax4,x,y,wpxy'/in.wpmax/w0,exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);
contour(ax5,x,z,wpxz'/in.wpmax/w0,exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);
contour(ax6,y,z,wpyz'/in.wpmax/w0,exp(-1)*[1, 1],'Color',[1 1 1],'linewidth',1.5);

end

if saveit,
    vidobj = VideoWriter(sprintf('%s%s_movie_scat.avi',datadir,runname));
    vidobj.FrameRate = 5;
    open(vidobj);
end

%% loop

for m = 1:in.numfiles,
    
    EYslicez = fread(Efile,[in.Ny in.Nx],'double');
    
    if size(EYslicez,1) < in.Ny,
        fclose('all');
        break;
    end
    
    EZslicez = fread(Efile,[in.Ny in.Nx],'double');
    EYslicey = fread(Efile,[in.Nz in.Nx],'double');
    EZslicey = fread(Efile,[in.Nz in.Nx],'double');
    EYslicex = fread(Efile,[in.Nz in.Ny],'double');
    EZslicex = fread(Efile,[in.Nz in.Ny],'double');
    Eyinc = fread(Efile,in.Nx,'double');
    Ezinc = fread(Efile,in.Nx,'double');

    JYslicez = fread(Jfile,[in.Ny in.Nx],'double');
    JZslicez = fread(Jfile,[in.Ny in.Nx],'double');
    JYslicey = fread(Jfile,[in.Nz in.Nx],'double');
    JZslicey = fread(Jfile,[in.Nz in.Nx],'double');
    JYslicex = fread(Jfile,[in.Nz in.Ny],'double');
    JZslicex = fread(Jfile,[in.Nz in.Ny],'double');
   
    HYslicez = fread(Hfile,[in.Ny in.Nx],'double');
    HZslicez = fread(Hfile,[in.Ny in.Nx],'double');
    HYslicey = fread(Hfile,[in.Nz in.Nx],'double');
    HZslicey = fread(Hfile,[in.Nz in.Nx],'double');
    HYslicex = fread(Hfile,[in.Nz in.Ny],'double');
    HZslicex = fread(Hfile,[in.Nz in.Ny],'double');
    Hyinc = fread(Hfile,in.Nx,'double');
    Hzinc = fread(Hfile,in.Nx,'double');
    
    Eyinc2D = repmat(Eyinc,1,in.Ny);
    Ezinc2D = repmat(Ezinc,1,in.Ny);
    % subtract incident field from inside to get pure scattered field
    
    EYscatz = EYslicez;
    EYscatz(by1:by2,bz1:bz2) = EYscatz(by1:by2,bz1:bz2) - Eyinc2D(by1:by2,bz1:bz2)';
    EYscaty = EYslicey;
    EYscaty(bx1:bx2,bz1:bz2) = EYscaty(bx1:bx2,bz1:bz2) - Eyinc2D(bx1:bx2,bz1:bz2)';
    EYscatx = EYslicex;
    EYscatx(bx1:bx2,by1:by2) = EYscatx(bx1:bx2,by1:by2) - Eyinc((in.Nx+1)/2 + 1);
    
    EYtotz = EYscatz + Eyinc2D';
    EYtoty = EYscaty + Eyinc2D';
    EYtotx = EYscatx;
    
    EZscatz = EZslicez;
    EZscatz(by1:by2,bz1:bz2) = EZscatz(by1:by2,bz1:bz2) - Ezinc2D(by1:by2,bz1:bz2)';
    EZscaty = EZslicey;
    EZscaty(bx1:bx2,bz1:bz2) = EZscaty(bx1:bx2,bz1:bz2) - Ezinc2D(bx1:bx2,bz1:bz2)';
    EZscatx = EZslicex;
    EZscatx(bx1:bx2,by1:by2) = EZscatx(bx1:bx2,by1:by2);
    
    EZtotz = EZscatz + Ezinc2D';
    EZtoty = EZscaty + Ezinc2D';
    EZtotx = EZscatx;
    
    set(im1,'CData',log10(abs(EYslicez)));
    set(im2,'CData',log10(abs(EYslicey)));
    set(im3,'CData',log10(abs(EYslicex)));
    set(ti1,'String',sprintf('Ey  (x-y slice) at time %.0f ns',(m-1)*writedt*1e9));
    set(ti2,'String',sprintf('Ey  (x-z slice) at time %.0f ns',(m-1)*writedt*1e9));
    set(ti3,'String',sprintf('Ey  (y-z slice) at time %.0f ns',(m-1)*writedt*1e9));
    
    if strcmp(plotfield2, 'total'),
    
        set(im4,'CData',log10(abs(EZslicez)));
        set(im5,'CData',log10(abs(EZslicey)));
        set(im6,'CData',log10(abs(EZslicex)));
        set(ti4,'String',sprintf('Ez  (x-y slice) at time %.0f ns',(m-1)*writedt*1e9));
        set(ti5,'String',sprintf('Ez  (x-z slice) at time %.0f ns',(m-1)*writedt*1e9));
        set(ti6,'String',sprintf('Ez  (y-z slice) at time %.0f ns',(m-1)*writedt*1e9));
    
    else
        
        set(im4,'CData',log10(abs(EYscatz)));
        set(im5,'CData',log10(abs(EYscaty)));
        set(im6,'CData',log10(abs(EYscatx)));
        set(ti4,'String',sprintf('Ey scattered field at time %.0f ns',(m-1)*writedt*1e9));
        set(ti5,'String',sprintf('Ey scattered field at time %.0f ns',(m-1)*writedt*1e9));
        set(ti6,'String',sprintf('Ey scattered field at time %.0f ns',(m-1)*writedt*1e9));
        
    end
    
    
    drawnow;
    if saveit,
        %print(h1,'-dpng',sprintf('%sframe-%03d.png',datadir,m));
        currframe = getframe(h1);
        writeVideo(vidobj,currframe);
    end
    pause(0.1);
    
end

if saveit,
    close(vidobj);
end

fclose('all');


