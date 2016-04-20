function glean_run(GLEAN)                                 
% Runs a GLEAN analysis.
%
% GLEAN = GLEAN_RUN(GLEAN)
%
% Adam Baker 2015


% Check settings and directories before running
GLEAN = glean_check(GLEAN);
save(GLEAN.name,'GLEAN')

pretty_string('RUNNING GLEAN ANALYSIS')

fprintf('Running GLEAN specified in: \n%s',GLEAN.name);

% Run the envelope state:
glean_timeseries(GLEAN)

% Run the subspace state:
glean_subspace(GLEAN)

% Run the model state:
glean_model(GLEAN)

pretty_string('GLEAN ANALYSIS COMPLETE')

end
