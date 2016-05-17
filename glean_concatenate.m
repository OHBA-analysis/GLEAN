function [dataConcat,subIndex,trlIndex,cndIndex,sessle] = glean_concatenate(GLEAN,freqOpt)
% Concatenates GLEAN subspace data into a single data matrix
%
% [dataConcat,subIndex] = GLEAN_CONCATENATE(GLEAN)
%
% REQUIRED INPUTS:
%   GLEAN       - An existing GLEAN analysis
%   freqOpt     - (optional) string to specify frequency band handling:
%                   - 'concatenate' concatenate frequencies
%                   - 'keep'        keep frequencies as third dimension
% OUTPUTS:
%   dataConcat  - [channels x timepoints (x frequency)] data matrix
%   subIndex    - [1 x timepoints] matrix of subject indicators
%
% Adam Baker 2015

% Concatenate subjects - should be the same as within hmm - so maybe put
% into function.

if nargin == 1
    freqOpt = 'keep';
end

trlIndex    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
cndIndex    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
subIndex    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
dataConcat  = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays

for session = 1:numel(GLEAN.data)
    D = spm_eeg_load(GLEAN.subspace.data{session});
    subIndex{session} = session * ones(1,D.nsamples);
    % Collate trialwise information
    if isfield(GLEAN.model.settings.hmm,'trialwise')
        
        %        trlIndex{session} = D.trl_indices + length(trlIndex{1});
        if session>1
            trlIndex{session} = D.trl_indices + trlIndex{session-1}(end);
        else
            trlIndex{session} = D.trl_indices ;
        end
        conds = reshape(D.cnd_indices,sum(D.trl_indices == 1),[]);
        cndIndex{session} = conds(1,:);
    end
    
    switch freqOpt
        case 'concatenate'
            dat = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
        case 'keep'
            dat = permute(D(:,:,:,:),[1 3 2 4]);
        otherwise
            error('Unknown freqOpt')
    end
    
    sessle{session}=nsamples(D)*ntrials(D);

    dataConcat{session} = dat;

end

if isfield(GLEAN.model.settings.hmm,'trialwise')
trlIndex    = cell2mat(trlIndex);
cndIndex    = cell2mat(cndIndex);
end

subIndex    = cell2mat(subIndex);
dataConcat  = cell2mat(dataConcat);

end
