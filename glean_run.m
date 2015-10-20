function glean_run(GLEAN)                                 
% Runs a GLEAN analysis.
%
% GLEAN = GLEAN_RUN(GLEAN)
%
% Adam Baker 2015


% Check settings and directories before running
GLEAN = glean_check(GLEAN);
save(GLEAN.name,'GLEAN')

% Run the envelope state:
glean_envelope(GLEAN)

% Run the subspace state:
glean_subspace(GLEAN)

% Run the model state:
glean_model(GLEAN)

% Run the results state:
glean_results(GLEAN)

end






                    
