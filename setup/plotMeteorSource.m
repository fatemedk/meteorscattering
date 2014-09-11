function [source,ax] = plotMeteorSource(in,ax)

loadconstants;

rdelayA = 1.6*sqrt(6)/pi/in.f0;
rdelayB = (1.6+1/sqrt(3)) * sqrt(6)/pi/in.f0;

t = 0:(in.tsteps-1);

source.A = zeros(size(t));
source.B = zeros(size(t));

if in.sourcetype == 1,
    
    source.A = sin(2*pi*in.f0*t*in.dt);
    source.B = cos(2*pi*in.f0*t*in.dt);

elseif in.sourcetype == 2,
    
    source.A = sin(2*pi*in.f0*(t-3*in.fsig)*in.dt) .* exp(-(t-3*in.fsig).^2/(2*in.fsig^2));
    source.B = cos(2*pi*in.f0*(t-3*in.fsig)*in.dt) .* exp(-(t-3*in.fsig).^2/(2*in.fsig^2));
    
elseif in.sourcetype == 3,
    
    source.A = (1 - 2*(pi*in.f0*(t*in.dt-rdelayA)).^2) * exp(-(pi*in.f0*(t*in.dt-rdelayA)).^2);
    source.B = (1 - 2*(pi*in.f0*(t*in.dt-rdelayB)).^2) * exp(-(pi*in.f0*(t*in.dt-rdelayB)).^2);

    
end


plot(ax,t*in.dt*1e9,source.A,'b');
hold(ax,'on');
plot(ax,t*in.dt*1e9,source.B,'r');