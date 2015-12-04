function res = glean_group_temporal_stats(GLEAN,ev,contrasts)
% Group differences in state temporal properties.
%
% GLEAN_GROUP_TEMPORAL_STATS(GLEAN,design,contrasts)%
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
%   ev        - Design matrix of regressors [sessions x regressors]
%   contrasts - Matrix of contrasts to compute [contrasts x regressors]
%
% OUTPUTS:
%   res       - results structure with the following fields for each
%               temporal statistic computed:
%                .label   - label for this statistic
%                .units   - units for this statistic
%                .stats   - [sessions x states] matrix of statistics
%                .pvalues - [contrasts x states] matrix of p-values
%
% Adam Baker 2015


hmmstats = glean_hmm_stats(GLEAN);

num_perms       = 1e4;
num_sessions    = size(ev,1);
num_states      = numel(hmmstats);
num_contrasts   = size(contrasts,1);

res.nOccurrences        = struct('label','Number of occurrences',...
                          'units', '');
res.FractionalOccupancy = struct('label','Fractional occupancy',...
                          'units', '%');
res.MeanLifeTime        = struct('label','Mean life time',...
                          'units', 'ms');
res.MeanIntervalLength  = struct('label','Mean interval length',...
                          'units', 'ms');
res.Entropy             = struct('label','Entropy',...
                          'units', '');

for stat = fieldnames(res)'
    
    % Compute t-stats for specified design matrix and contrasts:
    stats = cat(1,hmmstats.(char(stat)))';
    stats = group_imputation(stats,ev); % impute missing values:
    [~,~,~,tstats] = glean_glm(stats,ev,contrasts);
    
    % Permutation testing:
    permuted_tstats = zeros(num_perms,num_contrasts,num_states);
    for perm = 1:num_perms
        ev_perm = ev(randperm(num_sessions),:);
        stats = cat(1,hmmstats.(char(stat)))';
        stats = group_imputation(stats,ev_perm);
        [~,~,~,permuted_tstats(perm,:,:)] = glean_glm(stats,ev_perm,contrasts);
    end
    
    % p-values from permutations:
    pvalues = zeros(num_contrasts,num_states);
    for c = 1:num_contrasts
        for k = 1:num_states
            counts = sum(abs(tstats(c,k)) <= abs(permuted_tstats(:,c,k)));
            pvalues(c,k) = (counts + 1) / (num_perms + 1);
        end
    end
    res.(char(stat)).stats      = stats;
    res.(char(stat)).tstats     = tstats;
    res.(char(stat)).CI.lower   = squeeze(prctile(permuted_tstats,2.5)); 
    res.(char(stat)).CI.upper   = squeeze(prctile(permuted_tstats,97.5));  
    res.(char(stat)).pvalues    = pvalues;
    
end

end




