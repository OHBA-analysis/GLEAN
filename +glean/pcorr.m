function GLEAN = pcorr(GLEAN,settings)
% Compute partial correlation between state time courses and envelopes.
%
% glean.pcorr(GLEAN,settings)
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
%
% Adam Baker 2015


res = 'pcorr';

% Check inputs:
% ...


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

% Load model
model = load(GLEAN.model.model);

    
for subspace = cellstr(settings.space)
    
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
    
    % Regressors are the state time courses (HMM) or independent components (ICA)
    switch char(intersect(lower(fieldnames(GLEAN.model.settings)),{'hmm','ica'}));
        case 'hmm'
            regressors = cell2mat(arrayfun(@(k) model.hmm.statepath==k,1:model.hmm.K,'UniformOutput',0));
            session_maps = nan(num_channels,model.hmm.K,num_frequencies,num_sessions);
        case 'ica'
            regressors = model.ica.tICs';
            session_maps = nan(num_channels,size(model.ica.tICs,1),num_frequencies,num_sessions);
    end
    
    % Compute partial correlation map for each subject
    for session = 1:numel(GLEAN.data)
        
        if num_frequencies == 1
            session_maps(:,:,1,session) = glean.regress(GLEAN.(data).data{session},regressors(model.subIndx==session,:),res);
        else
            session_maps(:,:,:,session) = glean.regress(GLEAN.(data).data{session},regressors(model.subIndx==session,:),res);
        end
        % Save the session specific maps
        disp(['Saving ' char(subspace) ' partial correlation maps for session ' num2str(session)])
        
        for f = 1:num_frequencies
            map = session_maps(:,:,f,session);
            switch settings.format
                case 'mat'
                    save2mat(results.(char(subspace)).sessionmaps{session}{f},map);
                case 'nii'
                    save2nii(map,results.(char(subspace)).sessionmaps{session}{f},char(subspace))
            end
        end
        
    end
    % Compute group map as the average of the session maps
    group_maps = nanmean(session_maps,4);
    
    % Save the group averaged maps
    disp(['Saving ' char(subspace) ' group partial correlation map'])
    for f = 1:num_frequencies
        map = group_maps(:,:,f);
        switch settings.format
            case 'mat'
                save2mat(results.(char(subspace)).groupmaps{f},map);
            case 'nii'
                save2nii(map,results.(char(subspace)).groupmaps{f},char(subspace))
        end
    end
end


% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings  = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');



    function save2nii(map,fname,space)
        % Save a matrix MAP to a nifti file with filename FNAME using a mask
        % appropriate for the SPACE the map is in (voxelwise or parcelwise)
        switch space
            case 'voxel'
                glean.writenii(map,fname,GLEAN.envelope.settings.mask);
            case 'parcel'
                map = glean.parcellation2map(map,GLEAN.subspace.settings.parcellation.file,GLEAN.envelope.settings.mask);
                glean.writenii(map,fname,GLEAN.envelope.settings.mask);
        end
    end

    function save2mat(fname, map)                                          %#ok<INUSD>
        % save MAP as a .mat file
        dirName = fileparts(fname);
        if ~isdir(dirName),
            mkdir(dirName);
        end%if
        save(fname, 'map');
    end%save2mat

end











