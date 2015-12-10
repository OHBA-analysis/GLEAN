function GLEAN = glean_group_netmats(GLEAN,settings)
% Group differences in state network matrices.
%
% GLEAN_GROUP_NETMATS(GLEAN,settings)
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
%                                network matrix
%                 .pvalues     - [contrasts x states x channels x channels] 
%                                p-values
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
[dat,subIndx] = glean_concatenate(GLEAN);
dat = normalise(dat');
 
% Compute netmats and group differences for each state
for k = 1:num_states
    
    sprintf('Computing network matrix for state %i of %i \n',k,num_states)
    
    % Select data for this state:
    datk = dat(model.hmm.statepath == k,:);
    subk = subIndx(model.hmm.statepath == k);
            
    % Compute netmat and group difference
    [C,netmat] = compute_netmat_diff(settings.design);
    results.netmats{k} = netmat;
    
    % Permutation testing:
    C_permuted = zeros([num_perms,size(C)]);
    for perm = 1:num_perms
        permuted_design = settings.design(randperm(num_sessions),:);
        C_permuted(perm,:,:,:) = compute_netmat_diff(permuted_design);
    end
    
    % p-values from permutations:
    pvalues = zeros(num_contrasts,num_channels,num_channels);
    for con = 1:num_contrasts
        for i = 1:num_channels
            for j = 1:num_channels
                counts = sum(abs(C(con,i,j)) <= abs(C_permuted(:,con,i,j)));
                pvalues(con,i,j) = (counts + 1) / (num_perms + 1);
            end
        end
    end
   
    results.pvalues{k} = pvalues;
    
end

% Reshape outputs
results.netmats = permute(cat(4,results.netmats{:}),[1 4 2 3]);
results.pvalues = permute(cat(4,results.pvalues{:}),[1 4 2 3]);

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
