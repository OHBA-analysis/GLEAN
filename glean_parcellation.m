function D = glean_parcellation(S)
% Computes node time series from a MEEG object using a parcellation.
%
% D = GLEAN_PARCELLATION(S)
%
% REQUIRED INPUTS:
%
%   S.D                 - SPM12 MEEG object
%   S.parcellation      - .nii or .nii.gz file at the same gridstep as S.D
%                         OR a .mat file containing a matrix of dimensions
%                         voxels x parcels with the same number of voxels
%                         as S.D
%   S.mask              - mask for this parcellation
%   S.orthogonalisation - Orthogonalisation protocol to apply:
%                         ['none','symmetric','closest','householder']
%   S.method            - method for reconstructing parcel time course:
%                         ['PCA','mean','peakVoxel','spatialBasis']
%   S.trialwise         - indicate whether to compute and othogonalise the
%                         parcellation trial-by-trial or all at once [0 1]
% OUTPUTS:
%   D                   - Newly created SPM12 MEEG object
%
% Adam Baker 2015

D = spm_eeg_load(S.D);

S.orthogonalisation = ft_getopt(S,'orthogonalisation','none');
S.method            = ft_getopt(S,'method','PCA');
S.trialwise         = ft_getopt(S,'trialwise',0);

% Parcellation - pass in as a P*V matrix or 1*V vector
if ~isempty(strfind(S.parcellation,'.nii'))
    parcellation = readnii(S.parcellation,S.mask);
elseif ~isempty(strfind(S.parcellation,'.mat'))
    parcellation = load(S.parcellation);
    if length(fieldnames(parcellation)) == 1
        parcellation = parcellation.(char(fieldnames(parcellation)));
    else
        error('.mat file should contain only one variable')
    end
else
    error('S.parcellation needs to be a .nii or .mat file');
end

if S.trialwise == 1
    % Compute parcellation trial-by-trial
    nodedata = zeros(size(parcellation,2),D.nsamples,D.ntrials);
    for idx = 1:D.ntrials
        nodedata(:,:,idx) = ROInets.get_node_tcs(D(:,:,idx), parcellation, S.method);
        nodedata(:,:,idx) = ROInets.remove_source_leakage(nodedata(:,:,idx), S.orthogonalisation);
    end
else
    % Concatenate epochs and parcel/orthogonalise in one go
    dat = reshape(D(:,:,:),D.nchannels,D.nsamples*D.ntrials);
    nodedata = ROInets.get_node_tcs(dat, parcellation, S.method);
    clearvars dat; % this could be huge
    nodedata = ROInets.remove_source_leakage(nodedata, S.orthogonalisation);
    % restack into trials
    nodedata = reshape(nodedata,size(parcellation,2),D.nsamples,D.ntrials);
end

% Save data to new MEEG object
outfile = fullfile(D.path,D.fname);
Dnode = clone(montage(D,'switch',0),outfile,[size(nodedata,1),D.nsamples,D.ntrials]);
Dnode = chantype(Dnode,1:Dnode.nchannels,'VE');
Dnode(:,:,:) = nodedata;
D = Dnode; % For output
D.save;

end
