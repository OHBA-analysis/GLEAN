function P = state_transitions(GLEAN)
% Computes the state transition probabilities for each session
%
% P = glean.state_transitions(GLEAN)
%
% REQUIRED INPUTS:
%   GLEAN   - An existing GLEAN analysis
%
% OUTPUTS:
%   P       - cell array of [K*K] state transition probabilities
%
% Adam Baker 2015

model = load(GLEAN.model.model);

% Could replace this with the correct way, which is reevaluating the
% transition matrix probabilistically based on each subject's data -
% e.g. getTransProbMats in HMMMAR, but need to figure out what to do
% with group-level normalisation.

num_sessions = numel(GLEAN.data);

P = cell(1,num_sessions);
for session = 1:num_sessions
    
    sp_sub = model.hmm.statepath(model.subIndx == session);
    transitions = [sp_sub(1:end-1) sp_sub(2:end)];
    
    p = zeros(model.hmm.K);
    for i = 1:model.hmm.K
        for j = 1:model.hmm.K
            p(i,j) = sum(ismember(transitions,[i,j],'rows')) / ...
                sum(transitions(:,1) == i);
        end
    end
    
    P{session} = p;
    
end


