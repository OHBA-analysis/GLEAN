function [dataConcat,subIndex] = glean_concatenate(GLEAN)
% Concatenates GLEAN subspace data into a single data matrix
%
% [dataConcat,subIndex] = GLEAN_CONCATENATE(GLEAN)
%
% REQUIRED INPUTS:
%   GLEAN       - An existing GLEAN analysis
%
% OUTPUTS:
%   dataConcat  - [channels x timepoints] data matrix
%   subIndex    - [1 x timepoints] matrix of subject indicators
%
% Adam Baker 2015

% Concatenate subjects - should be the same as within hmm - so maybe put
% into function.

subIndex    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
dataConcat  = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays

for session = 1:numel(GLEAN.data)
    D = spm_eeg_load(GLEAN.subspace.data{session});
    subIndex{session} = session * ones(1,D.nsamples);
    dat = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
    dataConcat{session} = dat;
end

subIndex    = cell2mat(subIndex);
dataConcat  = cell2mat(dataConcat);

end