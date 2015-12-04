function GLEAN = glean_group_power(GLEAN,design,contrasts)
% Group differences in state-specific activity.
%
% GLEAN_GROUP_POWER(GLEAN,design,contrasts)
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   design    - Design matrix of regressors [sessions x regressors]
%   contrasts - Matrix of contrasts to compute [contrasts x regressors]
%
% Adam Baker 2015

res = 'group_power';

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

% Append settings to GLEAN:
GLEAN.results.settings.(res).format     = 'nii';
GLEAN.results.settings.(res).space      = 'voxel'; % until I figure out a way to make weights normalisation work with orthogonalisation
GLEAN.results.settings.(res).design     = design;
GLEAN.results.settings.(res).contrasts  = contrasts;

% Set up data paths:
GLEAN = glean_data(GLEAN);
save(GLEAN.name,'GLEAN');

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
    num_channels    = D.nchannels;
    num_frequencies = D.nfrequencies;
    
    % Compute session level mean of envelope
    M = zeros(num_sessions,num_channels,num_frequencies);
    for session = 1:num_sessions
        
        % Get weights normalisation:
        D = spm_eeg_load(GLEAN.data{session});
        montage1 = montage(D,'getmontage',1);
        montage2 = montage(D,'getmontage',2);
        weights_norm = montage2.tra(:,1)./montage1.tra(:,1);
        
        % Load envelope data
        D = spm_eeg_load(GLEAN.(data).data{session});
        [~,M(session,:)] = glean_cov(D);
        
        % Remove weights normalisation
        M(session,:) = M(session,:) ./ weights_norm';
        
        % Save session means
        for f = 1:num_frequencies
            niifile = GLEAN.results.(res).(char(subspace)).sessionmaps{session}{f};
            map = squeeze(M(session,:,f))';
            switch char(subspace)
                case 'voxel' % write as 4D
                    writenii(map,niifile,GLEAN.envelope.settings.mask);
                case 'parcel' % write as 2D
                    map = permute(map,[1,3,4,2]);
                    writenii(map,niifile);
            end
            disp(['Saving ' char(subspace) 'wise mean maps for session ' num2str(session)]);
        end
    end
    
    % Run FSL randomise to perform permutation testing & FWE correction
    for f = 1:num_frequencies

        % Save means
            input_nii = fullfile(tmpdir,sprintf('cope_%i.nii.gz',f));
            map = M(:,:,f)';
            switch char(subspace)
                case 'voxel' % write as 4D
                    writenii(map,input_nii,GLEAN.envelope.settings.mask);
                case 'parcel' % write as 2D
                    map = permute(map,[1,3,4,2]);
                    writenii(map,input_nii);
            end
            
            % Run randomise
            output_nii = fullfile(tmpdir,sprintf('randomise_%i',f));
            command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
                              input_nii, output_nii, design_file, contrast_file);
            [~,~] = system(command);
            
            switch char(subspace)
                case 'voxel' % read as 4D
                     FWE_corrected_tstats = readnii([output_nii,'_tstat1.nii'], ...
                                                 GLEAN.envelope.settings.mask);      
%                     FWE_corrected_tstats = readnii([output_nii,'_vox_corrp_tstat1.nii'], ...
%                                                 GLEAN.envelope.settings.mask);                
                case 'parcel' % read as 2D
                     FWE_corrected_tstats = readnii([output_nii,'_tstat1.nii']);                      
%                     FWE_corrected_tstats = readnii([output_nii,'_vox_corrp_tstat1.nii']);  
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











