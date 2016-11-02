function [map,var_e,dof] = regress(D,regressors,mode,normalisation)
% Create spatial maps via multiple regression of HMM or ICA time courses.
%
% D = glean.regress(D,regressors,mode)
%
% REQUIRED INPUTS:
%   D               - Name of an SPM12 MEEG object
%   regressors      - [samples x regressors] matrix of temporal regressors
%   mode            - Type of map to create:
%                       'cope'  - contrast of parameter estimates
%                       'tstat' - t-statistic
%                       'pcorr' - partial correlation
%   normalisation   - Optional normalisation factor for each voxel/channel 
%                     (e.g. to perform weights normalisation) 
% 
% OUTPUTS:
%   map         - [voxels x regressors (x frequency)] spatial map
%   var_e       - (optional) variance of the residuals
%   dof         - (optional) number of degrees of freedom

% Adam Baker 2015


D = spm_eeg_load(D);

F = D.nfrequencies;
if isempty(F)
    F = 1;
end

if nargin < 4 || strcmp(mode,'pcorr')
    normalisation = ones(D.nchannels,1);
end


% Remove empty regressors:
regressors = double(regressors);
reg2use = ~all(regressors==0);
regressors = regressors(:,reg2use);


nRegressors = size(regressors,2);

map = nan(D.nchannels,length(reg2use),F);

if strcmp(mode,'pcorr')
    % Make std = 1
    x = bsxfun(@rdivide,regressors,std(regressors));
else
    x = regressors;
end

if all(isinf(x(:))) % Maybe warn instead
    x(:) = 0;
end

pinvxtx = pinv(x'*x);
dof     = D.nsamples - nRegressors;

% It's slow(ish) to read data from D channel by channel so read in as much
% as possible:
blks = memblocks(size(D),1);

for iblk = 1:size(blks,1)
    
    for f = 1:F
        
        if strcmp(D.transformtype,'TF')
            Dblk = permute(D(blks(iblk,1):blks(iblk,2),f,:,1),[1 3 2]);
        else
            Dblk = D(blks(iblk,1):blks(iblk,2),:,1);
        end
        
        for v = blks(iblk,1):blks(iblk,2)
            
            y = Dblk(blks(iblk,1):blks(iblk,2) == v,:,1)';
            
            % Apply normalisation
            if strcmp(mode,'pcorr')
                y = (y - mean(y))./std(y);
            else
                y = y .* normalisation(blks(iblk,1):blks(iblk,2) == v);
                y = y - mean(y);
            end
            
            beta = pinvxtx * x' * y;                
            e = y - x*beta;
            var_e = diag(e'*e / dof);
            
            switch mode
                case {'pcorr','cope'}
                    map(v,reg2use,f) = beta;
                case 'tstat'
                    t = beta ./ diag(sqrt(var_e*pinvxtx));
                    map(v,reg2use,f) = t;
            end
            
        end
        
    end
    
end
