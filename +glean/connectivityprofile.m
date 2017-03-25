function GLEAN = connectivityprofile(GLEAN,settings)
% Computes the "connectivity profile" P for each HMM state. 
% The connectivity profile is defined for each state as the Euclidean 
% distance between a particular voxel's covariance with all other voxels
% within the state versus outside of the state.
%
% GLEAN = glean.connectivityprofile(GLEAN,SETTINGS)
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   settings  - Structure with the following fields:
%                 .format    - Format to save maps as
%                                permitted: 'mat','nii'
%                                default: 'mat'
%
% Adam Baker 2015

res = 'connectivity_profile';

if numel(GLEAN.envelope.settings.freqbands) > 1
    error('Connectivity profile is not yet supported for multiband data');
end

% Remove existing results:
if isfield(GLEAN.results,res)
    GLEAN.results = rmfield(GLEAN.results,res);
end

results = struct;

% Set up file names
results_dir = fullfile(GLEAN.results.dir,res); % make this an option

group_maps   = fullfile(results_dir,strcat('group_',res));

% Duplicate maps across each frequency band:
group_maps = strcat(group_maps,settings.format);
results.groupmaps = group_maps;


% Load model
model = load(GLEAN.model.model);




Nchannels = size(model.hmm.state(1).Cov,1);

C = zeros(Nchannels,Nchannels,model.hmm.K);

for k = 1:model.hmm.K
    C(:,:,k) = corrcov(model.hmm.state(k).Cov);
end

P = zeros(Nchannels,model.hmm.K);

channels = 1:Nchannels;
states = 1:model.hmm.K;

C(repmat(logical(eye(Nchannels)),[1,1,model.hmm.K])) = nan;

for channel = channels
    
    for state = states
        
        centroid = mean(C(channel,channels~=channel,states~=state),3);
        cp = C(channel,channels~=channel,state);
        
        P(channel,state) = norm((centroid - cp),2);
        
    end
    
end


% Save the group maps
map = P;
switch settings.format
    case 'mat'
        save(results.groupmaps,'map');
    case 'nii'
        map = glean.parcellation2map(map,GLEAN.subspace.settings.parcellation.file,GLEAN.envelope.settings.mask);
        glean.writenii(map,results.groupmaps,GLEAN.envelope.settings.mask);
end


% Append results and settings to GLEAN:
GLEAN.results.(res) = results;
GLEAN.results.(res).settings  = settings;

% Save updated GLEAN:
save(GLEAN.name,'GLEAN');

