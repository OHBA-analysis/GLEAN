function [outfile] = stacktrialsmeeg( D, outpath, conds );
%% function [outfile] = stacktrialsmeeg( D, outpath, conds );
%
% Function for concatenating epoched spm files into one continuous trial
%
% Inputs
% D: epoched spm object
% outpath: FULL path to output file
% cond: optional cell array of strings which could be used in a call to D.indtrial
%
% Outputs
% outfile: new, concatenated spm object


if nargin < 3 || isempty(conds)
    % Use all the trials
    condinds = 1:size(D,3);
else
    % Use only trials on specified condition
    condinds = D.indtrial(conds);
end

% Extract relevent condition labels
conditions = D.conditions(condinds);

%% Extract and reshape dataset
data = D(:,:,condinds);
dim_stats = size(data);
if ndims(dim_stats) < 4;dim_stats = [dim_stats 1]; end

% cont_data now contains one trial with all epoched stacked
cont_data = reshape(data,dim_stats(1),dim_stats(2)*dim_stats(3),dim_stats(4));
cont_dim_stats = [size(cont_data)];
if ndims(cont_dim_stats) < 3;cont_dim_stats = [cont_dim_stats 1]; end

% This could be pretty big...
clearvars data;

% The actual copy
clone(D, outpath, [ cont_dim_stats(1), cont_dim_stats(2), cont_dim_stats(3) ]);

% Create the trial and condition indices within the stacked dataset
trl_indices = zeros(1,cont_dim_stats(2)); % Indexes trial number
cnd_indices = zeros(1,cont_dim_stats(2)); % Indexes condition within D.condlist

count = 1;
for idx = 1:dim_stats(2):cont_dim_stats(2)
    cnd_indices(idx:idx+dim_stats(2)) = find(ismember(D.condlist',conditions{count}));
    trl_indices(idx:idx+dim_stats(2)) = count;
    count = count + 1;
end

% Save it out
outfile = spm_eeg_load(outpath);
outfile.trl_indices = trl_indices;
outfile.cnd_indices = cnd_indices;
outfile.save;
