function [beta, cope, varcope, tstat, dof] = glm(y,X,contrasts)
% Compute parameter estimates and contrasts for the GLM: Y = X * Beta + e
%
% [beta, cope, varcope, tstat] = glean.glm(y,X,contrasts)
%
% REQUIRED INPUTS:
%   y           - data matrix [observations x channels]
%   X           - design matrix of regressors [observations x regressors]
%   contrasts   - matrix of contrasts to compute [contrasts x regressors]
%
% OUTPUTS:
%   beta        - parameter estimates [regressors x channels]
%   cope        - contrasts of parameter estimates [contrasts x channels]
%   varcope     - variance for each contrasts [contrasts x channels]
%   tstat       - t-statistic for each contrast [contrasts x channels]
%   dof         - degrees of freedom
% Adam Baker 2015

% Check dimensions are consistent:
assert(size(y,1) == size(X,1),'y and X must have the same number of observations')

% Dimensions:
num_observations = size(X,1);
num_regressors   = size(X,2);
num_channels     = size(y,2);
if nargin > 2
    num_contrasts = numel(contrasts);
    assert(size(contrasts,2) == size(X,2),'contrasts and X must have the same number of regressors')
else
    num_contrasts = 0;
end

% Precompute some things:
pinvX   = pinv(X);
pinvXTX = pinv(X'*X);
dof     = num_observations - num_regressors;

% Preallocate:
beta    = zeros(num_regressors, num_channels);
var_e   = zeros(num_channels,1);

% Loop over blocks of channels
blks = osl_memblocks(size(y),2);
for i = 1:size(blks,1)
    
    blk_inds = blks(i,1):blks(i,2);
    
    % Parameter estimates:
    beta(:,blk_inds) = pinvX * y(:,blk_inds);
    
    % Error variance:
    e = y(:,blk_inds) - X * beta(:,blk_inds);
    var_e = diag(e'*e / dof);
 
end


% Compute contrasts:
if num_contrasts > 0
    
    % Contrast of parameter estimates:
    if nargout > 1
        cope = contrasts * beta;
    end
    
    % Variance of contrasts:
    if nargout > 2
        varcope = diag(contrasts * pinvXTX * contrasts') * var_e';
    end
    
    % t-statistics:
    if nargout > 3
        tstat = cope ./ sqrt(varcope);
    end
    
end





