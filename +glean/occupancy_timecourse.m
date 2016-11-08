function [occupancy,subIndx] = occupancy_timecourse(GLEAN,winsize,type)
% Compute fractional occupancy time courses from HMM state time courses
%
% [occupancy,subIndx] = glean.occupancy_timecourse(GLEAN,settings)
%
% Computes fractional occupancy time courses for each state, by applying a
% moving average window to the HMM state time courses
%
% REQUIRED INPUTS:
%   GLEAN     - An existing GLEAN analysis
%   winsize   - size of moving average window in samples
%
%
% OPTIONAL INPUTS:
%   type      - type of state path to use ['viterbi','gamma']
%               default - 'viterbi'
%
% OUTPUTS:
%   occupancy - occupancy time course 
%   subIndx   - session indices for occupancy time course
%
% Adam Baker 2015

model = load(GLEAN.model.model);

if ~exist(type,'var')
    type = 'viterbi';
end

window = rectwin(winsize) ./ winsize;

occupancy = zeros(model.hmm.K,ceil(length(model.hmm.statepath)/winsize));

for k = 1:model.hmm.K
    
    switch lower(char(type))
        case 'viterbi'
            statepath = double(model.hmm.statepath == k);
        case 'gamma'
            statepath = model.hmm.train.Gamma(:,k);
    end
    
    occupancy(k,:) = downsample(conv(statepath,window,'same'),winsize);
    
end

subIndx = downsample(model.subIndx,winsize);

end


