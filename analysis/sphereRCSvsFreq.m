% conducting sphere RCS, to compare my results...

function rcs = sphereRCSvsFreq(f,msize)

% inputs are frequency array f, and meteor size msize in meters

iterations = 100;
rcs = zeros(size(f));

c = 3e8;

for k = 1:length(f),
    
  freq = f(k);
  lambda = c / freq;     % calculate wavelength in meters for each frequency step.
  ka = 2 * pi * msize / lambda ;        % a frequently used constant {dimensionless}.
  s = sqrt( pi / 2 / ka ) ;         % to convert cartesian Bessel and Hankel functions to
                                    %   spherical Bessel and Hankel functions, increase the
                                    %   ORDER of the function by 1/2, and multiply by 's'.
  n = 1:iterations;                 % "n" is the ORDER of the Bessel and Hankel functions.

  B1 = besselj( n + 1/2, ka );
  B2 = besselh( n + 1/2, 2, ka );
  B3 = besselj( n + 1/2 - 1, ka );
  B4 = besselh( n + 1/2 - 1, 2, ka );

  a = ( s*B1 ) ./ ( s*B2 );
  b = ( ka * s*B3 - n .* s.*B1 ) ./ ( ka * s*B4 - n .* s.*B2 );
  rcs(k) = (lambda^2 / pi) * ( abs( sum( (-1).^n .* ( n + 1/2 ) .* ( b(n) - a(n) ) ) ) )^2;

end