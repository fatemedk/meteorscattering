% calculate electron-neutral, electron-electron, and other collision
% frequencies based on neutral density (altitude) and electron density

function in = meteorCollisions(in)

loadconstants;


%% start with electron neutral collisions - directly related to altitude

% uses as inputs in.altitude (in meters) and in.wp (array)
n0 = MSISatmosphere1(0);
nd = MSISatmosphere1(in.altitude/1e3);

% mobility is related to nd; then collision frequency to mobility

mue = 1.4856 * n0.total / nd.total;

in.coll.en = QE/(ME * mue);


%% okay, next. we want coulomb collisions (e-e, e-i). 

ne = in.wp.^2 * ME * e0 / QE^2;
ni = ne;

% need to assume temperatures. For now, let's assume it's cold, so Tn = Te
% = Ti taken from MSIS above. Note that colder = higher collision
% frequency!

in.Te = nd.temp * kB/QE;        % converted to eV
in.Ti = in.Te * kB/QE;

% Coulomb Logarithm: approximate by 20
lnLam = 20;

% e-e collision frequency from equation 2.86 of J. D. Callen, Fundamentals
% of Plasma Physics, 2006. 

% we only need electron-ion and electron-electron collisions

z = 1;      % assume singly-charged
in.coll.ei = 6.6e-11 * ni * z^2 * (MI/ME)^(-1/2) * in.Te^(-3/2) * (lnLam / 17);

% electron-electron collisions

in.coll.ee = 6.6e-11 * ni * z^2 * (ME/ME)^(-1/2) * in.Te^(-3/2) * (lnLam / 17);