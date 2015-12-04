function glean_group_state_lifetimes(GLEAN,design,contrasts)

res = 'group_state_lifetimes';

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

% Append settings to GLEAN:
GLEAN.results.settings.(res).format     = 'eps';
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

% glean_hmm_stats returns a 1 x K struct containing temporal statistics for
% each session
hmmstats = glean_hmm_stats(GLEAN);


stats = cat(1,hmmstats.MeanLifeTime);

% Impute missing values:


input_nii = fullfile(tmpdir,'tmp.nii.gz');
writenii(permute(stats,[1,3,4,2]),input_nii);

output_nii = fullfile(tmpdir,'randomise');
command = sprintf('randomise -i %s -o %s -d %s -t %s -x', ...
    input_nii, output_nii, design_file, contrast_file);
fprintf('Running FSL randomise')
[~,~] = system(command);


% Run FSL randomise to perform permutation testing (over subjects0 & 
% FWE correction (over states)
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
                FWE_corrected_tstats(:,k) = readnii([output_nii,'_vox_corrp_tstat1.nii'], ...
                    GLEAN.envelope.settings.mask);
            case 'parcel' % read as 2D
                FWE_corrected_tstats(:,k) = readnii([output_nii,'_vox_corrp_tstat1.nii']);
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












