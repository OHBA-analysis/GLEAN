function run(GLEAN)                                 
% Runs a GLEAN analysis.
%
% GLEAN = glean.run(GLEAN)
%
% Adam Baker 2015


% Check settings and directories before running
GLEAN = glean.check(GLEAN);
save(GLEAN.name,'GLEAN')

pretty_string('RUNNING GLEAN ANALYSIS')

fprintf('Running GLEAN specified in: \n%s',GLEAN.name);

% Run the envelope state:
glean.envelope(GLEAN)

% Run the subspace state:
glean.subspace(GLEAN)

% Run the model state:
glean.model(GLEAN)

pretty_string('GLEAN ANALYSIS COMPLETE')

end
