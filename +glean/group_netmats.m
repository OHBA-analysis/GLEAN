function GLEAN = group_netmats(GLEAN,settings)
% Group differences in state network matrices.
%
% glean.group_netmats(GLEAN,settings)
%
% Computes network matrices within each HMM state, by selecting only data
% that is within the state and group of interest. Group differences are
% computed between the networks matrices, and significance testing is
% achieved by random permutation of the group labels.
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   settings  - Structure with the following fields:
%                 .design    - [sessions x regressors] design matrix 
%                              (must contain only 0s and 1s)
%                 .contrasts - [contrasts x regressors] matrix of contrasts 
%                              to compute 
%
% OUTPUTS:
%   GLEAN     - An existing GLEAN analysis with new results field
%               GLEAN.results.group_transitions, containing the
%               following fields:
%                 .netmats     - [groups x states x channels x channels] 
%                                network matrix for each frequency band
%                 .pvalues     - [contrasts x states x channels x channels] 
%                                p-values for each frequency band
%
% Adam Baker 2015

res = 'group_netmats';

% Network matrix function
netfun = @(x) fisher(corr(x));

% Check inputs:
% ...

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

model = load(GLEAN.model.model);

results = struct;

num_perms       = 1e4;
num_sessions    = size(settings.design,1);
num_states      = model.hmm.K;
num_groups      = size(settings.design,2);
num_contrasts   = size(settings.contrasts,1);
num_channels    = size(model.hmm.state(1).Cov,1);
    
% Concatenate subjects
[dat,subIndx] = glean.concatenate(GLEAN);
dat = normalise(dat,2);
 
% Compute netmats and group differences for each state
for f = 1:size(dat,3)
    
    for k = 1:num_states
        
        sprintf('Computing network matrix for state %i of %i \n',k,num_states)
        
        % Select data for this state:
        datk = dat(:,model.hmm.statepath == k,f)';
        subk = subIndx(model.hmm.statepath == k);
        
        % Compute netmat and group difference
        [C,netmat] = compute_netmat_diff(settings.design);
        results.netmats{f}{k} = netmat;
        
        % Permutation testing:
        counts = zeros(num_contrasts,num_channels,num_channels);
        for perm = 1:num_perms
            permuted_design = settings.design(randperm(num_sessions),:);
            C_permuted = compute_netmat_diff(permuted_design);       
            counts = counts + double(abs(C) <= abs(C_permuted));
        end
        
        % p-values from permutations:
        pvalues = (counts + 1) ./ (num_perms + 1);
        
        results.pvalues{f}{k} = pvalues;
        
    end
    
end

% Reshape outputs
for f = 1:size(dat,3)
    results.netmats{f} = permute(cat(4,results.netmats{f}{:}),[1 4 2 3]);
    results.pvalues{f} = permute(cat(4,results.pvalues{f}{:}),[1 4 2 3]);
end
% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');



function [con,c] = compute_netmat_diff(design)
% Compute netmats for each group
sessions = 1:num_sessions;

c = zeros(num_groups,num_channels,num_channels);
for group = 1:size(design,2)
    % Get indices for this group, for this state
    group_inds = ismember(subk,sessions(design(:,group) == 1));
    datk_group = datk(group_inds,:);
    c(group,:,:) = netfun(datk_group);
end

% Compute contrasts
con = zeros(num_contrasts,num_channels,num_channels);
for ch = 1:num_channels
    con(:,:,ch) = settings.contrasts * c(:,:,ch);
end

end

end
