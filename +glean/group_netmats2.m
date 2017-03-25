function GLEAN = group_netmats2(GLEAN,settings)
% Group differences in state network matrices.
%
% glean.group_netmats2(GLEAN,settings)
%
% Computes network matrices within each HMM state for each subject by 
% selecting only data that is within the state and group of interest. 
% Group differences are computed between the networks matrices, and 
% significance testing is achieved by random permutation of the group 
% labels.
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
%                 .tstats      - [contrasts x states x channels x channels] 
%                                t-statistics for each frequency
%
% Adam Baker 2015

res = 'group_netmats';

% Network matrix function
netfun = @(x) corr(x);

% Check inputs:
% ...

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

model = load(GLEAN.model.model);

results = struct;

num_sessions    = size(settings.design,1);
num_states      = model.hmm.K;
num_channels    = size(model.hmm.state(1).Cov,1);
num_contrasts   = size(settings.contrasts,1);

% Create temporary directory
tmpdir = tempname;
mkdir(tmpdir);
c = onCleanup(@() system(['rm -r ' tmpdir]));

% Compute netmats and group differences for each state
for k = 1:num_states
    
    D = spm_eeg_load(GLEAN.subspace.data{1});
    
    C   = zeros(num_channels,num_channels,D.nfrequencies,num_sessions);
    dof = zeros(num_sessions,1);
    
    for session = 1:num_sessions
        D = spm_eeg_load(GLEAN.subspace.data{session}); 
        sp_sub = model.hmm.statepath(model.subIndx == session);
        dof(session) = sum(sp_sub == k);
        for f = 1:D.nfrequencies
            datk = permute(D(:,f,find(sp_sub == k)),[3,1,2]); %#ok
            C(:,:,f,session) = netfun(datk);
        end
    end
    
    
    % Remove sessions with dof less than nchannels 
    sessions2exclude = find(dof <= num_channels);
    warning('Excluding %i sessions from state %i due to too few samples (less than number of channels)', ...
            length(sessions2exclude),k);
    C(:,:,:,sessions2exclude) = [];
    design = settings.design;
    design(sessions2exclude,:) = [];
    disp(['State ' num2str(k) ' now has group sizes: ' num2str(sum(design))])
    
    
    % Save design matrix and contrast files
    design_file = fullfile(tmpdir,'design.mat');
    save_vest(design,design_file);
    contrast_file = fullfile(tmpdir,'design.con');
    save_vest(settings.contrasts,contrast_file);
    
    % Run FSL randomise to perform permutation testing & FWE correction
    randomise_output = zeros(num_channels,num_channels,D.nfrequencies,num_contrasts,num_states);
    % Save netmats for each state
    input_nii = fullfile(tmpdir,sprintf('netmat_%i.nii.gz',k));
    glean.writenii(C,input_nii);
     
    % Run randomise
    output_nii = fullfile(tmpdir,sprintf('randomise_%i',k));
    command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
                      input_nii, output_nii, design_file, contrast_file);
    fprintf('Running FSL randomise for state %i of %i\n',k,num_states)
    [~,~] = system(command);
        
    for con = 1:num_contrasts
        randomise_output(:,:,:,con,k) = glean.readnii([output_nii,'_tstat' num2str(con) '.nii']);
    end
    
end

% Reshape outputs
for f = 1:D.nfrequencies
    results.tstats{f} = permute(randomise_output(:,:,f,:,:),[4 5 1 2 3]);
end

% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');

end
