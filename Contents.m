%         ____ _     _____    _    _   _ 
%        / ___| |   | ____|  / \  | \ | |
%       | |  _| |   |  _|   / _ \ |  \| |
%       | |_| | |___| |___ / ___ \| |\  |
%        \____|_____|_____/_/   \_|_| \_|
%
%   Group Level Exploratory Analysis of Networks                              
%
%                Adam Baker 2015
%
%
% For details on how to set up a new GLEAN analysis type "help glean_setup"
% For details on particular GLEAN settings type "help glean_check"
% To run a GLEAN analysis type "run_glean(GLEAN)"
%
% If using this toolbox please consider citing the following:
%
% "Fast transient networks in spontaneous human brain activity", Baker el al., eLife, 2014
% "Spectrally resolved fast transient brain states in electrophysiological data", Vidaurre et al., NeuroImage, 2016
% "A symmetric multivariate leakage correction for MEG connectomes", Colclough et al., NeuroImage, 2015
%
% Primary functions:
%   glean.setup                   - Sets up a GLEAN analysis with particular filename, settings and data.
%   glean.check                   - Checks the settings and set up directory and data structures in GLEAN.
%   glean.directories             - Sets up the directory structure for the GLEAN analysis.
%   glean.data                    - Set up the directory structure for a new or existing GLEAN analysis.
%   glean.run                     - Runs a GLEAN analysis.
%   glean.envelope                - Runs the envelope stage of GLEAN.
%   glean.subspace                - Runs the subspace stage of GLEAN.
%   glean.model                   - Runs the model stage of GLEAN.
%   glean.path                    - Returns the path to the GLEAN code.
%   glean.load                    - Load an existing GLEAN analysis
%   glean.convert2spm             - Converts data to SPM12 format, saving the new .mat and .dat files.
%
% Main helper functions:
%   glean.hilbenv                 - Optimised Hilbert envelope computation for of MEEG data. 
%   glean.normalise               - Performs normalisation of MEEG data.
%   glean.concatenate             - Concatenates GLEAN subspace data into a single data matrix
%   glean.cov                     - Computes the covariance of a [channels x samples] efficiently.
%   glean.groupcov                - Efficiently compute a group covariance from multiple SPM files.
%   glean.parcellation            - Computes node time series from a MEEG object using a parcellation.
%   glean.regress                 - Create spatial maps via multiple regression of HMM or ICA time courses.
%   glean.infer_hmm               - Infers an hidden Markov model (HMM) with particular options.
%   glean.glm                     - Compute parameter estimates and contrasts for the GLM: Y = X * Beta + e
%   glean.variance                - Computes the temporal variance of data in MEEG object
%
% Outputs and statistics:
%   glean.hmm_stats               - Compute summary statistics from HMM statepath. 
%   glean.pcorr                   - Compute partial correlation between state time courses and envelopes.
%   glean.state_transitions       - Computes the state transition probabilities for each session
%   glean.connectivityprofile     - Computes the "connectivity profile" P for each HMM state. 
%   glean.occupancy_timecourse    - Compute fractional occupancy time courses from HMM state time courses
%   glean.temporal_stats_plot     - Pretty plotting of HMM statistics:
%   glean.temporal_stats          - State temporal properties.
%   glean.plot_timecourse         - Plots GLEAN inferred time courses.
%   glean.write_model_timecourses - Write GLEAN model time courses as SPM12 MEEG objects.
%
% Group level analyses:
%   glean.group_power             - Group differences in oscillatory power.
%   glean.group_lifetimes         - Group differences in state life times, across all states.
%   glean.group_netmats           - Group differences in state network matrices.
%   glean.group_netmats2          - Group differences in state network matrices, version 2.
%   glean.group_transitions       - Group differences in state transition probabilities.
%   glean.group_temporal_stats    - Group differences in state temporal properties.
%   glean.group_statewise_power   - Group differences in state-specific activity.
%   glean.group_state_distance    - Group differences in state covariance similarity.


