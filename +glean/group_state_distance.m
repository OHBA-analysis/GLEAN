function GLEAN = group_state_distance(GLEAN,settings)
% Group differences in state covariance similarity.
%
% GLEAN = GLEAN_GROUP_STATE_DISTANCE(GLEAN,settings)
%
% Computes group differences between the within-session state covariance
% similarity (based on the Riemannian distance between each state's
% covariance matrix). Group differences are given as p-values, computed via
% random permutation of the group labels.
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   settings  - Structure with the following fields:
%                 .design    - [sessions x regressors] design matrix
%                 .contrasts - [1 x regressors] matrix of contrast 
%                              to compute 
%
% Adam Baker 2015

res = 'group_state_distance';

% Check inputs:
% ...

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

model = load(GLEAN.model.model);
results = struct;

D = spm_eeg_load(GLEAN.subspace.data{1});

num_perms       = 1e4;
num_sessions    = size(settings.design,1);
num_states      = model.hmm.K;
num_groups      = size(settings.design,2);
num_contrasts   = size(settings.contrasts,1);
num_channels    = size(model.hmm.state(1).Cov,1);
num_frequencies = D.nfrequencies;
 
% Compute netmats and group differences for each state
for f = 1:num_frequencies
                    
        % Concatenate subjects
        [dat,subIndx] = glean.concatenate(GLEAN);
        dat = normalise(dat(:,:,f)');
        
        % Compute netmat and group difference
        C = group_netmats(settings.design);
        
        % Permutation testing:
        counts = zeros(num_contrasts,num_channels,num_channels);
        for perm = 1:num_perms
            permuted_design = settings.design(randperm(num_sessions),:);
            C_permuted = group_netmats(permuted_design);       
            counts = counts + double(abs(C) <= abs(C_permuted));
            
            if mod(perm,100) == 0
                fprintf('Permutation %i of %i \n',perm,num_perms)
            end
        end
        
        % p-values from permutations:
        pvalues = (counts + 1) ./ (num_perms + 1);
        
        results.pvalues{f}{k} = pvalues;
        
    
end

% Reshape outputs
for f = 1:size(dat,3)
    results.pvalues{f} = permute(cat(4,results.pvalues{f}{:}),[1 4 2 3]);
end
% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');



function con = group_netmats(design)
% Compute netmats for each group
sessions = 1:num_sessions;

R = zeros(num_groups,1);

for group = 1:size(design,2)
    % Get indices for this group, for this state
    group_inds = ismember(subIndx,sessions(design(:,group) == 1))';
    
    R(group) = 0;
    for i = 1:num_states
        datk_group = dat(group_inds == 1 & model.hmm.statepath == i,:);
        Ci = cov(datk_group);
        for j = i:num_states
            datk_group = dat(group_inds == 1 & model.hmm.statepath == j,:);
            Cj = cov(datk_group);
            R(group) = R(group) + riemann_dist(Ci,Cj);
        end
    end
    R(group) = R(group) ./ (0.5 * num_states * (num_states - 1));
    
end

% Compute contrasts
con = settings.contrasts * R;

end

end
