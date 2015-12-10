function GLEAN = glean_group_transitions(GLEAN,settings)
% Group differences in state transition probabilities.
%
% GLEAN_GROUP_TRANSITIONS(GLEAN,settings)
%
% Computes temporal statistics from the HMM state time courses for each
% subject and state, and tests for group differences in these statistics
% using permutation testing. The two stages of this analysis are:
%
% 1) Computation of temporal statistics (number of occurences, fractional
%    occupancy, mean lifetime and mean interval length) for each state and
%    session.
%
% 2) Group-level t-tests computed using the specified group design matrix 
%    and contrasts. P-values are computed for these statistics by
%    permutation of the group design matrix.
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   settings  - Structure with the following fields:
%                 .design    - [sessions x regressors] design matrix
%                 .contrasts - [contrasts x regressors] matrix of contrasts 
%                              to compute 
%
% OUTPUTS:
%   GLEAN     - An existing GLEAN analysis with new results field
%               GLEAN.results.group_transitions, containing the
%               following fields:
%                 .transitions - [sessions x states x states] transition
%                                matrix
%                 .tstats      - [contrasts x states x states] matrix of 
%                                t-statistics
%                 .pvalues     - [contrasts x states x states] matrix of 
%                                p-values
%
% Adam Baker 2015

res = 'group_transitions';

% Check inputs:
% ...

model = load(GLEAN.model.model);

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

results = struct;

num_perms       = 1e4;
num_sessions    = size(settings.design,1);
num_states      = model.hmm.K;
num_contrasts   = size(settings.contrasts,1);

% Compute session-wise transition probabilities
P = glean_state_transitions(GLEAN);

results.transitions = permute(cat(3,P{:}),[3 1 2]);

% Compute t-stats for specified design matrix and contrasts:
Pvec = reshape(cat(3,P{:}),num_states^2,num_sessions)';
Pvec = group_imputation(Pvec,settings.design); % impute missing values
[~,~,~,tstats] = glean_glm(Pvec,settings.design,settings.contrasts);

results.tstats  = reshape(tstats,num_states,num_states);

% Permutation testing:
permuted_tstats = zeros(num_perms,num_contrasts,num_states^2);
for perm = 1:num_perms
    permuted_design = settings.design(randperm(num_sessions),:);
    Pvec = reshape(cat(3,P{:}),num_states^2,num_sessions)';
    Pvec = group_imputation(Pvec,permuted_design); % impute missing values
    [~,~,~,permuted_tstats(perm,:,:)] = glean_glm(Pvec,permuted_design,settings.contrasts);
end

% p-values from permutations:
pvalues = zeros(num_contrasts,num_states^2);
for c = 1:num_contrasts
    for t = 1:num_states^2
        counts = sum(abs(tstats(c,t)) <= abs(permuted_tstats(:,c,t)));
        pvalues(c,t) = (counts + 1) / (num_perms + 1);
    end
end

results.pvalues = reshape(pvalues,num_states,num_states);


% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');

end

