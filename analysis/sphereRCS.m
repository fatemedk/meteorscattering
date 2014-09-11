% conducting sphere RCS, to compare my results...

% calculate and plot RCS of a perfectly conducting sphere using 
% Eq.(2.XXXXX) and produce plots similar to Figures (2.XXXX) and (2.XXXXX) 
% Spherical Bessel functions are computed using series approximation 
% and recursion. 

clear all; %close all;

eps   = 0.00001;
index = 0;

% ka limits are [0.01 - 100], log-spaced, 1000 points

kvalues = logspace(-1,log10(50),1000);


for kr = kvalues,
    
   index = index + 1;
   sphere_rcs   = 0. + 0.*1i;
   f1    = 0. + 1.*1i;
   f2    = 1. + 0.*1i;
   m     = 1.;
   n     = 0.;
   q     = -1.;
   
   % initially set del to huge value
   del =100000+100000*1i;
   while(abs(del) > eps)
      q   = -q;
      n   = n + 1;
      m   = m + 2;
      del = (2.*n-1) * f2 / kr-f1;
      f1  = f2;
      f2  = del;
      del = q * m /(f2 * (kr * f1 - n * f2));
      sphere_rcs = sphere_rcs + del;
   end
   sphere_rcsdb(index) = 20. * log10(abs(sphere_rcs));
   rcs(index)   = abs(sphere_rcs);
end

figure;
semilogx(kvalues,rcs);
xlabel('Sphere circumference in wavelengths');
ylabel('Normalized sphere RCS (sig / pi r^2)');
grid;
axis tight;

figure;
loglog(kvalues,rcs);
xlabel('Sphere circumference in wavelengths');
ylabel('Normalized sphere RCS (sig / pi r^2) - dB');
grid;
axis tight;