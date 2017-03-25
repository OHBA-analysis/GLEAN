function GLEAN = group_statewise_power(GLEAN,settings)
% Group differences in state-specific activity.
%
% GLEAN_GROUP_STATEWISE_POWER(GLEAN,settings)
%
% Computes familywise error corrected t-statistics between groups on the 
% state-specific differences in band-limited amplitude between. The two
% stages of this analysis are:
%
% 1) Multiple regression of the state time courses on the (weights
%    unnormalised*) band-limited envelope, for each session, resulting in
%    subject level contrasts of parameter estimates (COPEs) for each state.
%
% 2) Group-level tests using the specified group design matrix and
%    contrasts. This is performed using FSL randomise, to control the
%    family-wise error rate across voxels or parcels.
%
%    * NOTE - weights normalisation is NOT correctly implemented for the
%      orthogonalised parcellation case.
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   settings  - Structure with the following fields:
%                 .format    - Format to save maps as
%                                permitted: 'mat','nii'
%                                default: 'mat'
%                 .space     - Subspace to save maps as
%                                permitted: 'parcel','voxel',{'parcel','voxel'}
%                                default: 'voxel'
%                 .design    - [sessions x regressors] design matrix
%                 .contrasts - [1 x regressors] matrix of contrast 
%                              to compute 
%                 .FWEC      - Use FSL's family-wise error correction [0/1]
%                                default: 1
%
% Adam Baker 2015


res = 'group_state_power';

% Check inputs:
% ...
try
    FWEC = settings.FWEC == 1;
catch
    FWEC = 1;
end

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

results = struct;

% Set up file names
results_dir = fullfile(GLEAN.results.dir,res); % make this an option
[~,session_names] = cellfun(@fileparts,GLEAN.data,'UniformOutput',0);
for space = cellstr(settings.space)
    session_maps = fullfile(results_dir,char(space),strcat(session_names,'_',res));
    group_maps   = fullfile(results_dir,char(space),strcat('group_',res));
    
    % Duplicate maps across each frequency band:
    fstr = cellfun(@(s) regexprep(num2str(s),'\s+','-'), GLEAN.envelope.settings.freqbands,'UniformOutput',0);
    group_maps = strcat(group_maps,'_',fstr,'Hz.',settings.format);
    if ~isempty(session_maps)
        session_maps = cellfun(@(s) strcat(s,'_',fstr,'Hz.',settings.format),session_maps,'UniformOutput',0);
    end
    
    results.(char(space)).sessionmaps  = session_maps;
    results.(char(space)).groupmaps    = group_maps;
end

% Load HMM
model = load(GLEAN.model.model);

% Create temporary directory
tmpdir = tempname;
mkdir(tmpdir);
c = onCleanup(@() system(['rm -r ' tmpdir]));

% Save design matrix and contrast files
design_file = fullfile(tmpdir,'design.mat');
save_vest(settings.design,design_file);
contrast_file = fullfile(tmpdir,'design.con');
save_vest(settings.contrasts,contrast_file);
    
    
for subspace = cellstr(settings.space)
    
    switch char(subspace)
        case 'voxel'
            data = 'envelope';
        case 'parcel'
            data = 'subspace';
            % Warn if using orthogonalisation
            if ~strcmp(GLEAN.subspace.settings.parcellation.orthogonalisation,'none')
                warning('Weights normalisation as used by this function does not correctly account for the effects of orthogonalisation. This feature needs to be added.')
            end
    end
    
    % Load first session to get dimensionality
    D = spm_eeg_load(GLEAN.(data).data{1});
    
    num_sessions    = numel(GLEAN.data);
    num_states      = model.hmm.K;
    num_channels    = D.nchannels;
    num_frequencies = D.nfrequencies;
    num_contrasts   = size(settings.contrasts,1);
    
    % Compute session level COPEs using state time course
    copes = zeros(num_sessions,num_channels,num_states);
    for session = 1:num_sessions
 
        sp_sub = model.hmm.train.Gamma(model.subIndx == session,:);

        % Get weights normalisation:
        D = spm_eeg_load(GLEAN.data{session});
        montage1 = montage(D,'getmontage',1);
        montage2 = montage(D,'getmontage',2);
        weights_norm = nanmedian(montage2.tra ./ montage1.tra, 2);

        % Load envelope data
        D = spm_eeg_load(GLEAN.(data).data{session});
        
        % Weights un-normalised COPEs
        copes(session,:,:,:) = glean.regress(D,sp_sub,'cope',1./weights_norm);
        
        % Save session COPEs
        for f = 1:num_frequencies
            niifile = results.(char(subspace)).sessionmaps{session}{f};
            map = squeeze(copes(session,:,:,f));
            switch char(subspace)
                case 'voxel' % write as 4D
                    glean.writenii(map,niifile,GLEAN.envelope.settings.mask);
                case 'parcel' % write as 2D
                    map = permute(map,[1,3,4,2]);
                    glean.writenii(map,niifile);
            end
            disp(['Saving ' char(subspace) 'wise COPE maps for session ' num2str(session)]);
        end
    end
    
    % Run FSL randomise to perform permutation testing & FWE correction
    for f = 1:num_frequencies
        tstats = zeros(num_channels,num_states,num_contrasts);
        for k = 1:num_states
            % Save COPEs for each state
            input_nii = fullfile(tmpdir,sprintf('cope_%i_%i.nii.gz',k,f));
            map = copes(:,:,k,f)';
            switch char(subspace)
                case 'voxel' % write as 4D
                    glean.writenii(map,input_nii,GLEAN.envelope.settings.mask);
                case 'parcel' % write as 2D
                    map = permute(map,[1,3,4,2]);
                    glean.writenii(map,input_nii);
            end
            
            % Run randomise
            output_nii = fullfile(tmpdir,sprintf('randomise_%i_%i',k,f));
            command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
                              input_nii, output_nii, design_file, contrast_file);
            fprintf('Running FSL randomise for state %i of %i\n',k,model.hmm.K)
            [~,~] = system(command);
            
            switch char(subspace)
                case 'voxel' % read as 4D
                    if FWEC
                        tstats(:,k) = glean.readnii([output_nii,'_vox_corrp_tstat1.nii'],GLEAN.envelope.settings.mask);
                    else
                        tstats(:,k) = glean.readnii([output_nii,'_tstat1.nii'],GLEAN.envelope.settings.mask);
                    end
                case 'parcel' % read as 2D
                    if FWEC
                        tstats(:,k) = glean.readnii([output_nii,'_vox_corrp_tstat1.nii']);
                    else
                        tstats(:,k) = glean.readnii([output_nii,'_tstat1.nii']);
                    end
            end
        end
        
        
        
        
        
        % Write FWE corrected t-stats to group maps
        if strcmp(subspace,'parcel')
            tstats = glean.parcellation2map(tstats, GLEAN.subspace.settings.parcellation.file, GLEAN.envelope.settings.mask);
        end
        glean.writenii(tstats, results.(char(subspace)).groupmaps{f}, GLEAN.envelope.settings.mask);
    end
end


% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings  = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');

end











