function hmmstats = glean_hmm_stats(GLEAN)
% Compute summary statistics from HMM statepath. 
%
% HMMSTATS = GLEAN_HMM_STATS(GLEAN) returns a number of statistics from the
% HMM state path (fractional occupancy, number of occurences, mean life
% time, mean interval length, life times, interval lengths).
%
% Adam Baker 2015

model = load(GLEAN.model.model);
if ~isfield(model,'hmm')
    error([mfilename,' only supports HMM models'])
end

% Get sampling rate from enveloped data
method = char(intersect(fieldnames(GLEAN.subspace.settings),{'pca','parcellation','voxel'}));
if strcmp(method,'parcellation')
    D = spm_eeg_load(GLEAN.subspace.data{1});
else
    D = spm_eeg_load(GLEAN.envelope.data{1});
end
fsample = D.fsample;
clear D


hmmstats = struct;
for k = 1:model.hmm.K
    for s = 1:length(unique(model.subIndx))
        sp      = model.hmm.statepath(model.subIndx==s);
        gamma   = model.hmm.train.Gamma(model.subIndx==s,k);
        
        lifetimes = diff(logical2epoch(sp==k),[],2)./fsample;
        lifetimes = lifetimes(lifetimes~=0);
        
        intervals = diff(logical2epoch(sp~=k),[],2)./fsample;
        intervals = intervals(intervals~=0);
        
        hmmstats(k).nOccurrences(s)         = length(lifetimes);
        hmmstats(k).FractionalOccupancy(s)  = 100 * sum(sp==k) / length(sp); % percent
        hmmstats(k).MeanLifeTime(s)         = 1000*mean(lifetimes); % ms
        hmmstats(k).MeanIntervalLength(s)   = 1000*mean(intervals); % ms
        hmmstats(k).Entropy(s)              = entropy(gamma);
        hmmstats(k).LifeTimes{s}            = lifetimes;
        hmmstats(k).Intervals{s}            = intervals;
    end
end


end


function ev = logical2epoch(l,t)
% Convert a logical vector to [2 x N] array of [onsets, offsets]
l = l(:)';

if nargin < 2
    t = 1:length(l);
end

onsets  = t(diff([0 l]) ==  1);
offsets = t(diff([l 0]) == -1);

ev = [onsets; offsets]';

end



