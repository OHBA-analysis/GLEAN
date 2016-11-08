function stat = group_imputation(stat,groups)
% Columnwise imputation of missing values by assigning the group mean.
%
% stat = GROUP_IMPUTATION(stat,groups)
%
% REQUIRED INPUTS:
%   stat        - [subjects x channels] matrix of statistics to clean
%   groups      - [subjects x groups] matrix in one-hot format
%
% OUTPUTS:
%   stat        - [subjects x channels] matrix of statistics where missing
%                 values (NaNs) are replaced by the mean value within the
%                 relevant group.
%
% Adam Baker 2015

% Loop through all channels
for k = 1:size(stat,2)
    
    % Compute group means
    group_means = nansum(groups .* repmat(stat(:,k),1,size(groups,2))) ...
                         ./ sum(groups(~isnan(stat(:,k)),:));
    
    % Replace missing values with mean for the relevant group                  
    inds2impute = find(isnan(stat(:,k)));
    stat(inds2impute,k) = sum(groups(inds2impute,:) .* repmat(group_means,size(inds2impute)),2);
end

end
