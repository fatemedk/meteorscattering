% Set up random set of meteor simulations, with random size, shape,
% density.

clear all; close all;

loadconstants;

% set defaults
meteordefaults;

% anything you want to change?
inputs.submitjob = 1;

% master directory for set of runs
toprundir = '/shared/users/ram80/meteor/runs/randomset/';

% set of values to vary for meteor
sizerange = [0.1 1];        % meters
shaprange = [1 1];
densrange = [0.1 10];

numsims = 100;

% create vectors of random values
allsize = sizerange(1) + (sizerange(2)-sizerange(1))*rand(numsims,1);
allshap = shaprange(1) + (shaprange(2)-shaprange(1))*rand(numsims,1);
alldens = densrange(1) + (densrange(2)-densrange(1))*rand(numsims,1);

% run random simulations

for m = 1:numsims,
    
    % change variables as requested
    inputs.metsize = allsize(m);
    inputs.metshape = allshap(m);
    inputs.wpmax = alldens(m);
    
    inputs.runname = sprintf('random_%03d',m);
    inputs.rundir = [toprundir inputs.runname];
    
    % launch job
    [in,jobid] = meteorrun(inputs);
    
end
