% GLEAN DEMO SCRIPT FOR INFERRING A GROUP HMM FROM SOURCE SPACE MEG DATA
% Adam Baker, Jan 2016

% Directory of the data:
demo_dir = './data/GLEAN_demo/'; 

% List of data files. These data consist of 2 sessions of 10 minutes of 
% eyes open resting state data recorded from 6 subjects on an Elekta
% Neuromag system (using the 204 planar gradiometers only). The data were
% band-pass filtered between 4 and 30 Hz and projected onto a regular 8 mm
% grid in source space using an LCMV beamformer. The data are saved as
% SPM12 MEEG objects in sensor space, and the beamformer weights are saved
% as a virtual montage.
data = {'/data/s01_rest1.mat'
        '/data/s01_rest2.mat'
        '/data/s02_rest1.mat'
        '/data/s02_rest2.mat'
        '/data/s03_rest1.mat'
        '/data/s03_rest2.mat'
        '/data/s04_rest1.mat'
        '/data/s04_rest2.mat'
        '/data/s05_rest1.mat'
        '/data/s05_rest2.mat'
        '/data/s06_rest1.mat'
        '/data/s06_rest2.mat'};
data = fullfile(demo_dir,data);
    

% Directory of the GLEAN analysis:
glean_dir  = fullfile(demo_dir,'analysis'); 

% Name for this GLEAN analysis:
glean_name = fullfile(glean_dir,'glean_demo.mat');
    
% Clear settings
settings = struct;

% Envelope settings:
settings.envelope.overwrite = 0;
settings.envelope.log       = 0;
settings.envelope.fsample   = 20;
settings.envelope.mask      = fullfile(demo_dir,'/masks/MNI152_T1_8mm_brain.nii.gz');

% Subspace settings:
settings.subspace.overwrite                         = 0;
settings.subspace.normalisation                     = 'none';
settings.subspace.parcellation.file                 = fullfile(demo_dir,'/masks/fMRI_parcellation_ds8mm.nii.gz');
settings.subspace.parcellation.orthogonalisation    = 'symmetric';
settings.subspace.parcellation.method               = 'spatialBasis';

% Model settings:
settings.model.overwrite   = 0;
settings.model.hmm.nstates = 8;
settings.model.hmm.nreps   = 1;

% Set up the GLEAN settings, data paths etc:
GLEAN = glean_setup(glean_name,data,settings);

% Run the analysis:
glean_run(GLEAN)



