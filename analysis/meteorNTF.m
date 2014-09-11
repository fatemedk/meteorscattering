% processing of meteor NTF transformation

function sig = meteorNTF(datadir,doplot)

loadconstants;

load([datadir 'inputs.mat']);


%% load from meteor NTF output file

fid = fopen([datadir 'output_NTF.dat'],'r');
nntf = fread(fid,1,'int');  % number of frequencies
ntfb = fread(fid,1,'int');
fntf = fread(fid,nntf,'double');

xsize = in.Nx;
ysize = in.Ny;
zsize = in.Nz;
fbsize = ysize*zsize*2*nntf;
tbsize = xsize*zsize*2*nntf;
lrsize = xsize*ysize*2*nntf;

% DFTs of source function
Eydft = fread(fid,nntf*2,'double');
Ezdft = fread(fid,nntf*2,'double');

% equivalent currents

My.front = fread(fid,fbsize,'double');
Mz.front = fread(fid,fbsize,'double');
Jy.front = fread(fid,fbsize,'double');
Jz.front = fread(fid,fbsize,'double');
My.back = fread(fid,fbsize,'double');
Mz.back = fread(fid,fbsize,'double');
Jy.back = fread(fid,fbsize,'double');
Jz.back = fread(fid,fbsize,'double');

Mx.top = fread(fid,tbsize,'double');
Mz.top = fread(fid,tbsize,'double');
Jx.top = fread(fid,tbsize,'double');
Jz.top = fread(fid,tbsize,'double');
Mx.bot = fread(fid,tbsize,'double');
Mz.bot = fread(fid,tbsize,'double');
Jx.bot = fread(fid,tbsize,'double');
Jz.bot = fread(fid,tbsize,'double');

Mx.left = fread(fid,lrsize,'double');
My.left = fread(fid,lrsize,'double');
Jx.left = fread(fid,lrsize,'double');
Jy.left = fread(fid,lrsize,'double');
Mx.right = fread(fid,lrsize,'double');
My.right = fread(fid,lrsize,'double');
Jx.right = fread(fid,lrsize,'double');
Jy.right = fread(fid,lrsize,'double');

fclose(fid);

% 16 = ntfb+1;
% 186 = in.Nx-ntfb;
% 187 = in.Nx-ntfb+1;


%% now I need to reshape all of these arrays and fix 1/2's on edges

My.front = reshape(My.front,2*nntf,in.Nz,in.Ny);
% now just take area between ntf boundaries (16, 187 for now): z then y!
My.front = My.front(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
% corrections to edges! fields on edges are associated with 1/2 face area
My.front(:,:,1) = My.front(:,:,1) * 0.5;
My.front(:,:,end) = My.front(:,:,end) * 0.5;

My.back = reshape(My.back,2*nntf,in.Nz,in.Ny);
My.back = My.back(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
My.back(:,:,1) = My.back(:,:,1) * 0.5;
My.back(:,:,end) = My.back(:,:,end) * 0.5;

Jz.front = reshape(Jz.front,2*nntf,in.Nz,in.Ny);
Jz.front = Jz.front(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Jz.front(:,:,1) = Jz.front(:,:,1) * 0.5;
Jz.front(:,:,end) = Jz.front(:,:,end) * 0.5;

Jz.back = reshape(Jz.back,2*nntf,in.Nz,in.Ny);
Jz.back = Jz.back(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Jz.back(:,:,1) = Jz.back(:,:,1) * 0.5;
Jz.back(:,:,end) = Jz.back(:,:,end) * 0.5;

% next set of four fields: Mz and Jy on front and back

Mz.front = reshape(Mz.front,2*nntf,in.Nz,in.Ny);
Mz.front = Mz.front(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Mz.front(:,1,:) = Mz.front(:,1,:) * 0.5;
Mz.front(:,end,:) = Mz.front(:,end,:) * 0.5;

Mz.back = reshape(Mz.back,2*nntf,in.Nz,in.Ny);
Mz.back = Mz.back(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Mz.back(:,1,:) = Mz.back(:,1,:) * 0.5;
Mz.back(:,end,:) = Mz.back(:,end,:) * 0.5;

Jy.front = reshape(Jy.front,2*nntf,in.Nz,in.Ny);
Jy.front = Jy.front(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Jy.front(:,1,:) = Jy.front(:,1,:) * 0.5;
Jy.front(:,end,:) = Jy.front(:,end,:) * 0.5;

Jy.back = reshape(Jy.back,2*nntf,in.Nz,in.Ny);
Jy.back = Jy.back(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Jy.back(:,1,:) = Jy.back(:,1,:) * 0.5;
Jy.back(:,end,:) = Jy.back(:,end,:) * 0.5;

% next set of fields on top and bottom (y faces)

Mx.top = reshape(Mx.top,2*nntf,in.Nz,in.Nx);
Mx.top = Mx.top(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Mx.top(:,:,1) = Mx.top(:,:,1) * 0.5;
Mx.top(:,:,end) = Mx.top(:,:,end) * 0.5;

Mx.bot = reshape(Mx.bot,2*nntf,in.Nz,in.Ny);
Mx.bot = Mx.bot(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Mx.bot(:,:,1) = Mx.bot(:,:,1) * 0.5;
Mx.bot(:,:,end) = Mx.bot(:,:,end) * 0.5;

Jz.top = reshape(Jz.top,2*nntf,in.Nz,in.Nx);
Jz.top = Jz.top(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Jz.top(:,:,1) = Jz.top(:,:,1) * 0.5;
Jz.top(:,:,end) = Jz.top(:,:,end) * 0.5;

Jz.bot = reshape(Jz.bot,2*nntf,in.Nz,in.Ny);
Jz.bot = Jz.bot(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Jz.bot(:,:,1) = Jz.bot(:,:,1) * 0.5;
Jz.bot(:,:,end) = Jz.bot(:,:,end) * 0.5;

% next set of four: Mz and Jx on top and bottom

Mz.top = reshape(Mz.top,2*nntf,in.Nz,in.Ny);
Mz.top = Mz.top(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Mz.top(:,1,:) = Mz.top(:,1,:) * 0.5;
Mz.top(:,end,:) = Mz.top(:,end,:) * 0.5;

Mz.bot = reshape(Mz.bot,2*nntf,in.Nz,in.Ny);
Mz.bot = Mz.bot(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Mz.bot(:,1,:) = Mz.bot(:,1,:) * 0.5;
Mz.bot(:,end,:) = Mz.bot(:,end,:) * 0.5;

Jx.top = reshape(Jx.top,2*nntf,in.Nz,in.Ny);
Jx.top = Jx.top(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Jx.top(:,1,:) = Jx.top(:,1,:) * 0.5;
Jx.top(:,end,:) = Jx.top(:,end,:) * 0.5;

Jx.bot = reshape(Jx.bot,2*nntf,in.Nz,in.Ny);
Jx.bot = Jx.bot(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Jx.bot(:,1,:) = Jx.bot(:,1,:) * 0.5;
Jx.bot(:,end,:) = Jx.bot(:,end,:) * 0.5;

% next set of fields on left and right (z faces)

Mx.left = reshape(Mx.left,2*nntf,in.Ny,in.Nx);
Mx.left = Mx.left(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Mx.left(:,:,1) = Mx.left(:,:,1) * 0.5;
Mx.left(:,:,end) = Mx.left(:,:,end) * 0.5;

Mx.right = reshape(Mx.right,2*nntf,in.Ny,in.Nx);
Mx.right = Mx.right(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Mx.right(:,:,1) = Mx.right(:,:,1) * 0.5;
Mx.right(:,:,end) = Mx.right(:,:,end) * 0.5;

Jy.left = reshape(Jy.left,2*nntf,in.Ny,in.Nx);
Jy.left = Jy.left(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Jy.left(:,:,1) = Jy.left(:,:,1) * 0.5;
Jy.left(:,:,end) = Jy.left(:,:,end) * 0.5;

Jy.right = reshape(Jy.right,2*nntf,in.Ny,in.Nx);
Jy.right = Jy.right(:,(ntfb+1):(in.Nx-ntfb),(ntfb+1):(in.Nx-ntfb+1));
Jy.right(:,:,1) = Jy.right(:,:,1) * 0.5;
Jy.right(:,:,end) = Jy.right(:,:,end) * 0.5;

% last set of four: My and Jx on left and right

My.left = reshape(My.left,2*nntf,in.Ny,in.Nx);
My.left = My.left(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
My.left(:,1,:) = My.left(:,1,:) * 0.5;
My.left(:,end,:) = My.left(:,end,:) * 0.5;

My.right = reshape(My.right,2*nntf,in.Ny,in.Nx);
My.right = My.right(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
My.right(:,1,:) = My.right(:,1,:) * 0.5;
My.right(:,end,:) = My.right(:,end,:) * 0.5;

Jx.left = reshape(Jx.left,2*nntf,in.Ny,in.Nx);
Jx.left = Jx.left(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Jx.left(:,1,:) = Jx.left(:,1,:) * 0.5;
Jx.left(:,end,:) = Jx.left(:,end,:) * 0.5;

Jx.right = reshape(Jx.right,2*nntf,in.Ny,in.Nx);
Jx.right = Jx.right(:,(ntfb+1):(in.Nx-ntfb+1),(ntfb+1):(in.Nx-ntfb));
Jx.right(:,1,:) = Jx.right(:,1,:) * 0.5;
Jx.right(:,end,:) = Jx.right(:,end,:) * 0.5;


%% now combine into three complex structures, one for each frequency


Myf.front = 1i*My.front(1:2:end,:,:) + My.front(2:2:end,:,:);
Mzf.front = 1i*Mz.front(1:2:end,:,:) + Mz.front(2:2:end,:,:);
Jyf.front = 1i*Jy.front(1:2:end,:,:) + Jy.front(2:2:end,:,:);
Jzf.front = 1i*Jz.front(1:2:end,:,:) + Jz.front(2:2:end,:,:);
Myf.back = 1i*My.back(1:2:end,:,:) + My.back(2:2:end,:,:);
Mzf.back = 1i*Mz.back(1:2:end,:,:) + Mz.back(2:2:end,:,:);
Jyf.back = 1i*Jy.back(1:2:end,:,:) + Jy.back(2:2:end,:,:);
Jzf.back = 1i*Jz.back(1:2:end,:,:) + Jz.back(2:2:end,:,:);

Mxf.top = 1i*Mx.top(1:2:end,:,:) + Mx.top(2:2:end,:,:);
Mzf.top = 1i*Mz.top(1:2:end,:,:) + Mz.top(2:2:end,:,:);
Jxf.top = 1i*Jx.top(1:2:end,:,:) + Jx.top(2:2:end,:,:);
Jzf.top = 1i*Jz.top(1:2:end,:,:) + Jz.top(2:2:end,:,:);
Mxf.bot = 1i*Mx.bot(1:2:end,:,:) + Mx.bot(2:2:end,:,:);
Mzf.bot = 1i*Mz.bot(1:2:end,:,:) + Mz.bot(2:2:end,:,:);
Jxf.bot = 1i*Jx.bot(1:2:end,:,:) + Jx.bot(2:2:end,:,:);
Jzf.bot = 1i*Jz.bot(1:2:end,:,:) + Jz.bot(2:2:end,:,:);

Mxf.left = 1i*Mx.left(1:2:end,:,:) + Mx.left(2:2:end,:,:);
Myf.left = 1i*My.left(1:2:end,:,:) + My.left(2:2:end,:,:);
Jxf.left = 1i*Jx.left(1:2:end,:,:) + Jx.left(2:2:end,:,:);
Jyf.left = 1i*Jy.left(1:2:end,:,:) + Jy.left(2:2:end,:,:);
Mxf.right = 1i*Mx.right(1:2:end,:,:) + Mx.right(2:2:end,:,:);
Myf.right = 1i*My.right(1:2:end,:,:) + My.right(2:2:end,:,:);
Jxf.right = 1i*Jx.right(1:2:end,:,:) + Jx.right(2:2:end,:,:);
Jyf.right = 1i*Jy.right(1:2:end,:,:) + Jy.right(2:2:end,:,:);

Ezsource = 1i*Ezdft(1:2:end) + Ezdft(2:2:end);
Eysource = 1i*Eydft(1:2:end) + Eydft(2:2:end);


%% create 2D arrays (front, back, left, right, top, bottom) of distance from
% each grid point to observation point: just x-distance!

dist.front = ntfb*in.dx;
dist.back = (in.Nx-ntfb)*in.dx;
dist.top1 = repmat((ntfb:(in.Nx-ntfb))*in.dx,in.Nz-2*ntfb,1);
dist.bot1 = dist.top1;
dist.left1 = repmat((ntfb:(in.Nx-ntfb))*in.dx,in.Ny-2*ntfb,1);
dist.right1 = dist.left1;

dist.top2 = repmat((ntfb:(in.Nx-ntfb-1))*in.dx,in.Nz-2*ntfb+1,1) + in.dx/2;
dist.bot2 = dist.top2;
dist.left2 = repmat((ntfb:(in.Nx-ntfb-1))*in.dx,in.Ny-2*ntfb+1,1) + in.dx/2;
dist.right2 = dist.left2;


%% Okay. Now I just need to integrate over the faces with dA and exp(jkR).

Lx = zeros(nntf,1);
Ly = zeros(nntf,1);
Lz = zeros(nntf,1);
Nx = zeros(nntf,1);
Ny = zeros(nntf,1);
Nz = zeros(nntf,1);


for m = 1:nntf,
    k = 2*pi*fntf(m)/vp;
    
    % Lx, with contributions from Mx (top, bottom, left, right)
    temp1 = squeeze(Mxf.top(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.top1);
    temp2 = squeeze(Mxf.bot(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.bot1);
    temp3 = squeeze(Mxf.left(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.left1);
    temp4 = squeeze(Mxf.right(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.right1);
    Lx(m) = sum(temp1(:)) + sum(temp2(:)) + sum(temp3(:)) + sum(temp4(:));
    
    % Ly, with contributions from My (front, back, left, right)
    temp1 = squeeze(Myf.front(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.front);
    temp2 = squeeze(Myf.back(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.back);
    temp3 = squeeze(Myf.left(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.left2);
    temp4 = squeeze(Myf.right(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.right2);
    Ly(m) = sum(temp1(:)) + sum(temp2(:)) + sum(temp3(:)) + sum(temp4(:));
    
    % Lz, with contributions from Mz (front, back, top, bottom)
    temp1 = squeeze(Mzf.front(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.front);
    temp2 = squeeze(Mzf.back(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.back);
    temp3 = squeeze(Mzf.top(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.top2);
    temp4 = squeeze(Mzf.bot(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.bot2);
    Lz(m) = sum(temp1(:)) + sum(temp2(:)) + sum(temp3(:)) + sum(temp4(:));
    
    % Nx, with contributions from Jx (top, bottom, left, right)
    temp1 = squeeze(Jxf.top(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.top2);
    temp2 = squeeze(Jxf.bot(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.bot2);
    temp3 = squeeze(Jxf.left(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.left2);
    temp4 = squeeze(Jxf.right(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.right2);
    Nx(m) = sum(temp1(:)) + sum(temp2(:)) + sum(temp3(:)) + sum(temp4(:));
    
    % Ny, with contributions from Jy (front, back, left, right)
    temp1 = squeeze(Jyf.front(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.front);
    temp2 = squeeze(Jyf.back(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.back);
    temp3 = squeeze(Jyf.left(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.left1);
    temp4 = squeeze(Jyf.right(m,:,:)) * in.dx*in.dy .* exp(1j*k*dist.right1);
    Ny(m) = sum(temp1(:)) + sum(temp2(:)) + sum(temp3(:)) + sum(temp4(:));
    
    % Nz, with contributions from Jz (front, back, top, bottom)
    temp1 = squeeze(Jzf.front(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.front);
    temp2 = squeeze(Jzf.back(m,:,:)) * in.dy*in.dz .* exp(1j*k*dist.back);
    temp3 = squeeze(Jzf.top(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.top1);
    temp4 = squeeze(Jzf.bot(m,:,:)) * in.dx*in.dz .* exp(1j*k*dist.bot1);
    Nz(m) = sum(temp1(:)) + sum(temp2(:)) + sum(temp3(:)) + sum(temp4(:));
    
end


%% now convert to Ltheta, Lphi, Ntheta, Nphi

% this is trivial, because for our geometry, th = -z and ph = -y directions

Lth = -Lz;
Lph = -Ly;
Nth = -Nz;
Nph = -Ny;

Ethsource = -Ezsource;
Ephsource = -Eysource;

Pinc = (abs(Ethsource).^2 + abs(Ephsource).^2)/2/eta0;

k = 2*pi*fntf/vp;

%% finally, compute cross sections

% from lecture notes by Gedney
sig.thth = k.^2/(4*pi) .* abs(Lph + eta0*Nth).^2 ./ abs(Ethsource).^2;
sig.phth = k.^2/(4*pi) .* abs(Lth - eta0*Nph).^2 ./ abs(Ethsource).^2;
sig.phph = k.^2/(4*pi) .* abs(Lth - eta0*Nph).^2 ./ abs(Ephsource).^2;
sig.thph = k.^2/(4*pi) .* abs(Lph + eta0*Nth).^2 ./ abs(Ephsource).^2;

% equation 8.37 from Taflove, v2, page 366
sig.sig  = k.^2/(8*pi*eta0) .* (abs(Lph + eta0*Nth).^2 + abs(Lth - eta0*Nph).^2) ./ abs(Pinc);



%% compare to flat plate scatterer from theory

if doplot,
    
    h = figure;
    plot(fntf/1e6,20*log10(sig.sig),'bx-');
    hold on;
    plot(fntf/1e6,20*log10(sig.thth),'rx-');
    plot(fntf/1e6,20*log10(sig.phph),'ro-');
    
    if in.doplate == 1,
        sig.expected = (4*pi*in.platesize^4).*fntf.^2/vp^2;
        plot(fntf/1e6,20*log10(sig.expected),'kx-');
        outfile = sprintf('sig_plate_%02.0f',in.platesize);
        save(outfile,'in','sig');
    elseif in.doplate == 2,
        sig.expected = sphereRCSvsFreq(fntf,in.platesize);
        plot(fntf/1e6,20*log10(sig.expected),'kx-');
    end
    
    
    xlabel('Frequency (MHz)');
    ylabel('Radar cross section (dB m^2)');
    legend('NTF (by power)','NTF (th-th)','NTF (ph-ph)','Theory');
    
end

save([datadir 'rcs.mat'],'sig');
