function GLEAN = glean_group_statewise_power(GLEAN,design,contrasts)
% Group differences in state-specific activity.
%
% GLEAN_GROUP_STATEWISE_POWER(GLEAN,design,contrasts)
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
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   design    - Design matrix of regressors [sessions x regressors]
%   contrasts - Matrix of contrasts to compute [contrasts x regressors]
%
% Adam Baker 2015

res = 'group_state_power';

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

% Append settings to GLEAN:
GLEAN.results.settings.(res).format     = 'nii';
GLEAN.results.settings.(res).space      = GLEAN.results.settings.pcorr.space;
GLEAN.results.settings.(res).design     = design;
GLEAN.results.settings.(res).contrasts  = contrasts;

% Set up data paths:
GLEAN = glean_data(GLEAN);
save(GLEAN.name,'GLEAN');

% Load HMM
model = load(GLEAN.model.model);

% Create temporary directory
tmpdir = tempname;
mkdir(tmpdir);
c = onCleanup(@() system(['rm -r ' tmpdir]));

% Save design matrix and contrast files
design_file = fullfile(tmpdir,'design.mat');
save_vest(design,design_file);
contrast_file = fullfile(tmpdir,'design.con');
save_vest(contrasts,contrast_file);

% Warn if using orthogonalisation
if isfield(GLEAN.subspace.settings,'parcellation') && ~strcmp(GLEAN.subspace.settings.parcellation.orthogonalisation,'none')
    warning('Weights normalisation as used by this function does not correctly account for the effects of orthogonalisation. This feature needs to be added.')
end
    
    
for subspace = fieldnames(GLEAN.results.(res))'
    
    switch char(subspace)
        case 'voxel'
            data = 'envelope';
        case 'parcel'
            data = 'subspace';
    end
    
    % Load first session to get dimensionality
    D = spm_eeg_load(GLEAN.(data).data{1});
    
    num_sessions    = numel(GLEAN.data);
    num_states      = model.hmm.K;
    num_channels    = D.nchannels;
    num_frequencies = D.nfrequencies;
    
    % Compute session level COPEs using state time course
    copes = zeros(num_sessions,num_channels,num_states);
    for session = 1:num_sessions
        
        % Get session statepath:
%         sp_sub = cell2mat(arrayfun(@(k) ...
%                           model.hmm.statepath(model.subIndx == session) == k, ...
%                           1:num_states, ...
%                           'UniformOutput', 0));
        sp_sub = model.hmm.train.Gamma(model.hmm.statepath(model.subIndx == session),:);


        % Get weights normalisation:
        D = spm_eeg_load(GLEAN.data{session});
        montage1 = montage(D,'getmontage',1);
        montage2 = montage(D,'getmontage',2);
        weights_norm = montage2.tra(:,1)./montage1.tra(:,1);
        
        % Load envelope data
        D = spm_eeg_load(GLEAN.(data).data{session});
        
        % Weights un-normalised COPEs
        copes(session,:,:,:) = glean_regress(D,sp_sub,'cope',1./weights_norm);
        
        % Save session COPEs
        for f = 1:num_frequencies
            niifile = GLEAN.results.(res).(char(subspace)).sessionmaps{session}{f};
            map = squeeze(copes(session,:,:,f));
            switch char(subspace)
                case 'voxel' % write as 4D
                    writenii(map,niifile,GLEAN.envelope.settings.mask);
                case 'parcel' % write as 2D
                    map = permute(map,[1,3,4,2]);
                    writenii(map,niifile);
            end
            disp(['Saving ' char(subspace) 'wise COPE maps for session ' num2str(session)]);
        end
    end
    
    % Run FSL randomise to perform permutation testing & FWE correction
    for f = 1:num_frequencies
        FWE_corrected_tstats = zeros(num_channels,num_states);
        for k = 1:num_states
            % Save COPEs for each state
            input_nii = fullfile(tmpdir,sprintf('cope_%i_%i.nii.gz',k,f));
            map = copes(:,:,k,f)';
            switch char(subspace)
                case 'voxel' % write as 4D
                    writenii(map,input_nii,GLEAN.envelope.settings.mask);
                case 'parcel' % write as 2D
                    map = permute(map,[1,3,4,2]);
                    writenii(map,input_nii);
            end
            
            % Run randomise
            output_nii = fullfile(tmpdir,sprintf('randomise_%i_%i',k,f));
            command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
                              input_nii, output_nii, design_file, contrast_file);
            fprintf('Running FSL randomise for state %i of %i\n',k,model.hmm.K)
            [~,~] = system(command);
            
            switch char(subspace)
                case 'voxel' % read as 4D
                     FWE_corrected_tstats(:,k) = readnii([output_nii,'_tstat1.nii'], ...
                                                 GLEAN.envelope.settings.mask);      
%                     FWE_corrected_tstats(:,k) = readnii([output_nii,'_vox_corrp_tstat1.nii'], ...
%                                                 GLEAN.envelope.settings.mask);                
                case 'parcel' % read as 2D
                     FWE_corrected_tstats(:,k) = readnii([output_nii,'_tstat1.nii']);                      
%                     FWE_corrected_tstats(:,k) = readnii([output_nii,'_vox_corrp_tstat1.nii']);  
            end
        end
        
        % Write FWE corrected t-stats to group maps
        if strcmp(subspace,'parcel')
            FWE_corrected_tstats = parcellation2map(FWE_corrected_tstats, ...
                                                    GLEAN.subspace.settings.parcellation.file, ...
                                                    GLEAN.envelope.settings.mask);
        end
        writenii(FWE_corrected_tstats, ...
                 GLEAN.results.(res).(char(subspace)).groupmaps{f}, ...
                 GLEAN.envelope.settings.mask);
    end
end

end











