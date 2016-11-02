function GLEAN = group_temporal_stats(GLEAN,settings)
% Group differences in state temporal properties.
%
% GLEAN_GROUP_TEMPORAL_STATS(GLEAN,settings)
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
%               And optional fields:
%                 .plot      - [1/0] whether to create figures (default 0)
%                 .grouplbls - [sessions x 1] cell array of group labels
%                 .conlbls   - [contrasts x 1] cell array of contrast
%                              labels
%
% OUTPUTS:
%   GLEAN     - An existing GLEAN analysis with new results field
%               GLEAN.results.group_temporal_stats, containing the
%               following fields:
%                 .label   - label for this statistic
%                 .units   - units for this statistic
%                 .stats   - [sessions x states] matrix of statistics
%                 .pvalues - [contrasts x states] matrix of p-values
%
% Adam Baker 2015

res = 'group_temporal_stats';

% Check inputs:
% ...


% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

hmmstats = glean.hmm_stats(GLEAN);

num_perms       = 1e4;
num_sessions    = size(settings.design,1);
num_states      = numel(hmmstats);
num_contrasts   = size(settings.contrasts,1);

results = struct;
results.nOccurrences        = struct('label','Number of occurrences',...
                                     'units', '');
results.FractionalOccupancy = struct('label','Fractional occupancy',...
                                     'units', '%');
results.MeanLifeTime        = struct('label','Mean life time',...
                                     'units', 'ms');
results.MeanIntervalLength  = struct('label','Mean interval length',...
                                     'units', 'ms');
results.Entropy             = struct('label','Entropy',...
                                     'units', '');

for stat = fieldnames(results)'
    
    fprintf('Computing statistics for %s \n',results.(char(stat)).label)
    
    % Compute t-stats for specified design matrix and contrasts:
    stats = cat(1,hmmstats.(char(stat)))';
    stats = group_imputation(stats,settings.design); % impute missing values:
    [~,~,~,tstats] = glean.glm(stats,settings.design,settings.contrasts);
    
    % Permutation testing:
    permuted_tstats = zeros(num_perms,num_contrasts,num_states);
    for perm = 1:num_perms
        permuted_design = settings.design(randperm(num_sessions),:);
        stats = cat(1,hmmstats.(char(stat)))';
        stats = group_imputation(stats,permuted_design);
        [~,~,~,permuted_tstats(perm,:,:)] = glean.glm(stats,permuted_design,settings.contrasts);
    end
    
    % p-values from permutations:
    pvalues = zeros(num_contrasts,num_states);
    for c = 1:num_contrasts
        for k = 1:num_states
            counts = sum(abs(tstats(c,k)) <= abs(permuted_tstats(:,c,k)));
            pvalues(c,k) = (counts + 1) / (num_perms + 1);
        end
    end
    results.(char(stat)).stats      = stats;
    results.(char(stat)).tstats     = tstats;
    results.(char(stat)).CI.lower   = squeeze(prctile(permuted_tstats,2.5)); 
    results.(char(stat)).CI.upper   = squeeze(prctile(permuted_tstats,97.5));  
    results.(char(stat)).pvalues    = pvalues;
    
    if isfield(settings,'plot') && settings.plot == 1
        results_dir = fullfile(GLEAN.results.dir,res,char(stat));
        if ~isdir(results_dir)
            mkdir(results_dir);
        end
        results.(char(stat)).plots.stats = fullfile(results_dir,'groups.fig');
        results.(char(stat)).plots.tstats = fullfile(results_dir,'tstats.fig');
        
        glean.temporal_stats_plot(results.(char(stat)),settings);
    end
    
end


% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings  = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');

end




