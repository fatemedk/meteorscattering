% default inputs for emp 2D jobs

inputs.exefile = 'meteor';
inputs.exedir = '/shared/users/ram80/meteor/source/';

inputs.submitjob = 1;
inputs.cluster = 'batch';
inputs.numnodes = '8';

% do you want to save simulation movies?
inputs.savefields = 1;

% do you want to use the PML boundary?
inputs.dopml = 1;

% source frequency
inputs.f0 = 150e6;

% source gaussian width 
inputs.fsig = 200;

% source type: 
% 1 = sinusoid at f0;
% 2 = modulated gaussian with sigma fsig;
% 3 = ricker wavelet at f0
inputs.sourcetype = 2;

% source polarization: 1 = linear Ez, 2 = linear Ey, 3 = 45 degrees, 4 = RHCP, 5 = LHCP
inputs.sourcepol = 3;

% meteorsize, in meters
inputs.metsize = 0.5;
inputs.metshape = 1;        % relative size of two dimensions of meteor; 1 = circular, greater than 1 = elongated away from source

% number of grid cells per wavelength
inputs.gridfactor = 200;

% maximum plasma frequency of meteor (in wp/w0)
inputs.wpmax = 1;

% maximum collision frequency (NO LONGER USED)
inputs.numax = 1e7;

% altitude: used to determine neutral density and thus e-n collision freq.
inputs.altitude = 100e3;

% meteor location (0,0,0 is center of grid)
inputs.metloc = [0 0 0];

% meteor angle w.r.t source, in x-y plane, in degrees
inputs.metangle = 0;

% meteor type: gaussian = 0, other = 1 (not yet implemented)
inputs.mettype = 0;

% number of output files (time steps)
inputs.numfiles = 100;

% grid size
inputs.Nx = 201;
inputs.Ny = 201;
inputs.Nz = 201;

% plate scatterer 1 meter^2 in size
% if doplate = 1, does a plate scatterer; if doplate = 2, spherical conductor
inputs.doplate = 0;
inputs.platesize = 1;

% set up probe points. Define locations in km, then determine grid values

inputs.nprobes = 1;
inputs.probex = (inputs.Nx+1)/2;
inputs.probey = (inputs.Ny+1)/2;
inputs.probez = (inputs.Nz+1)/2;

% near to far frequencies

inputs.fntf = linspace(0.9,1.1,15) * inputs.f0; 
inputs.nntf = length(inputs.fntf);

% total time steps to run is the time to cover the field of view
% (Nx*dx/vp/dt) times this tfactor, should be greater than 1 obviously.

inputs.tfactor = 5;

% use magnetic field? Forces Lee & Kalluri method, much slower

inputs.doB0 = 0;

inputs.Bx = 0;
inputs.By = 0;
inputs.Bz = 5e-5;
