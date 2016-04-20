function [outfile] = stacktrialsmeeg( D, outpath, conds )
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
    condinds = 1:D.ntrials;
else
    % Use only trials on specified condition
    condinds = D.indtrial(conds);
end

% Extract relevent condition labels
conditions = D.conditions(condinds);

%% Extract and reshape dataset
if strcmp(D.transformtype,'TF')
    data = D(:,:,:,condinds);
    dim_stats = size(data);
elseif strcmp(D.transformtype,'time')
    data = D(:,:,condinds);
    dim_stats = [size(data) 1];
else
    error('MEEG transformtype not recognised');
end

%% cont_dim_stats contains the size of the dataset to output
if strcmp(D.transformtype,'TF')
    cont_data = reshape(data,dim_stats(1),dim_stats(2),dim_stats(3)*dim_stats(4));
    cont_dim_stats = size(cont_data);
elseif strcmp(D.transformtype,'time');
    cont_data = reshape(data,dim_stats(1),dim_stats(2)*dim_stats(3));
    cont_dim_stats = [size(cont_data,1) 1 size(cont_data,2)];
end
if numel(cont_dim_stats) == 3;cont_dim_stats = [cont_dim_stats 1]; end

% This could be pretty big...
clearvars data;

%% Create new dataset
clone(D, outpath, cont_dim_stats);
outfile = spm_eeg_load(outpath);

if strcmp(D.transformtype,'TF')
    outfile(:,:,:,:) = cont_data;
elseif strcmp(D.transformtype,'time')
    outfile(:,:,:) = permute(cont_data,[1,3,2]);
end

%% Create the trial and condition indices within the stacked dataset
trl_indices = zeros(1,cont_dim_stats(2)); % Indexes trial number
cnd_indices = zeros(1,cont_dim_stats(2)); % Indexes condition within D.condlist

count = 1;
for idx = 1:D.nsamples:outfile.nsamples
    cnd_indices(idx:idx+D.nsamples-1) = find(ismember(D.condlist',conditions{count}));
    trl_indices(idx:idx+D.nsamples-1) = count;
    count = count + 1;
end

%% Save it out
outfile.trl_indices = trl_indices;
outfile.cnd_indices = cnd_indices;
outfile.save;
