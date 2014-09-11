%% launch batch of EMP 2D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
meteordefaults;

% anything you want to change?
inputs.submitjob = 1;
inputs.savefields = 1;
inputs.gridfactor = 100;
inputs.doB0 = 1;

% master directory for set of runs
toprundir = '/shared/users/ram80/meteor/runs/magtest/';

% variable for batch of runs. name must match an input!
var1.name = 'metsize';
var1.values = [1];


% submit jobs

for m = 1:length(var1.values),
    
    % change variables as requested
    evalstr = ['inputs.' var1.name ' = ' num2str(var1.values(m)) ';'];
    eval(evalstr);
    
    inputs.runname = [var1.name '_' sprintf('%.2g',var1.values(m))];
    inputs.runname(strfind(inputs.runname,'+')) = '';
    inputs.rundir = [toprundir inputs.runname];
    
    % launch job
    [in,jobid] = meteorrun(inputs);
    
end
