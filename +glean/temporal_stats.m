function GLEAN = temporal_stats(GLEAN,settings)
% State temporal properties.
%
% glean.temporal_stats(GLEAN,settings)
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   settings  - Structure with the following fields:
%                 .plot      - [1/0] whether to create figures (default 0)
%
% OUTPUTS:
%   GLEAN     - An existing GLEAN analysis with new results field
%               GLEAN.results.group_temporal_stats, containing the
%               following fields:
%                 .label   - label for this statistic
%                 .units   - units for this statistic
%                 .stats   - [sessions x states] matrix of statistics
%
% Adam Baker 2015

res = 'temporal_stats';

% Check inputs:
% ...


% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

hmmstats = glean.hmm_stats(GLEAN);

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
    
    results.(char(stat)).stats = cat(1,hmmstats.(char(stat)))'; 
    
    if isfield(settings,'plot') && settings.plot == 1
        results_dir = fullfile(GLEAN.results.dir,res,char(stat));
        if ~isdir(results_dir)
            mkdir(results_dir);
        end
        results.(char(stat)).plots.stats = fullfile(results_dir,'stats.fig');       
        glean.temporal_stats_plot(results.(char(stat)),settings);
    end
    
end


% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings  = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');

end




