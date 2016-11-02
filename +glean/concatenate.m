function [dataConcat,subIndex] = concatenate(GLEAN,freqOpt)
% Concatenates GLEAN subspace data into a single data matrix
%
% [dataConcat,subIndex] = glean.concatenate(GLEAN)
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
    
subIndex    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
dataConcat  = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays

for session = 1:numel(GLEAN.data)
    D = spm_eeg_load(GLEAN.subspace.data{session});
    subIndex{session} = session * ones(1,D.nsamples);
    
    switch freqOpt
        case 'concatenate'
            dat = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
        case 'keep'
            dat = permute(D(:,:,:,:),[1 3 2 4]);
        otherwise
            error('Unknown freqOpt')
    end
    dataConcat{session} = dat;

end

subIndex    = cell2mat(subIndex);
dataConcat  = cell2mat(dataConcat);

end