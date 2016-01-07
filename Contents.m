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
%   glean_setup                   - Sets up a GLEAN analysis with particular filename, settings and data.
%   glean_check                   - Checks the settings and set up directory and data structures in GLEAN.
%   glean_directories             - Sets up the directory structure for the GLEAN analysis.
%   glean_data                    - Set up the directory structure for a new or existing GLEAN analysis.
%   glean_run                     - Runs a GLEAN analysis.
%   glean_envelope                - Runs the envelope stage of GLEAN.
%   glean_subspace                - Runs the subspace stage of GLEAN.
%   glean_model                   - Runs the model stage of GLEAN.
%   glean_path                    - Returns the path to the GLEAN code.
%   glean_load                    - Load an existing GLEAN analysis
%   glean_convert2spm             - Converts data to SPM12 format, saving the new .mat and .dat files.
%
% Main helper functions:
%   glean_hilbenv                 - Optimised Hilbert envelope computation for of MEEG data. 
%   glean_normalise               - Performs normalisation of MEEG data.
%   glean_concatenate             - Concatenates GLEAN subspace data into a single data matrix
%   glean_cov                     - Computes the covariance of a [channels x samples] efficiently.
%   glean_groupcov                - Efficiently compute a group covariance from multiple SPM files.
%   glean_parcellation            - Computes node time series from a MEEG object using a parcellation.
%   glean_regress                 - Create spatial maps via multiple regression of HMM or ICA time courses.
%   glean_infer_hmm               - Infers an hidden Markov model (HMM) with particular options.
%   glean_glm                     - Compute parameter estimates and contrasts for the GLM: Y = X * Beta + e
%   glean_variance                - Computes the temporal variance of data in MEEG object
%
% Outputs and statistics:
%   glean_hmm_stats               - Compute summary statistics from HMM statepath. 
%   glean_pcorr                   - Compute partial correlation between state time courses and envelopes.
%   glean_state_transitions       - Computes the state transition probabilities for each session
%   glean_connectivityprofile     - Computes the "connectivity profile" P for each HMM state. 
%   glean_occupancy_timecourse    - Compute fractional occupancy time courses from HMM state time courses
%   glean_temporal_stats_plot     - Pretty plotting of HMM statistics:
%   glean_temporal_stats          - State temporal properties.
%   glean_plot_timecourse         - Plots GLEAN inferred time courses.
%   glean_write_model_timecourses - Write GLEAN model time courses as SPM12 MEEG objects.
%
% Group level analyses:
%   glean_group_power             - Group differences in oscillatory power.
%   glean_group_lifetimes         - Group differences in state life times, across all states.
%   glean_group_netmats           - Group differences in state network matrices.
%   glean_group_netmats2          - Group differences in state network matrices, version 2.
%   glean_group_transitions       - Group differences in state transition probabilities.
%   glean_group_temporal_stats    - Group differences in state temporal properties.
%   glean_group_statewise_power   - Group differences in state-specific activity.
%   glean_group_state_distance    - Group differences in state covariance similarity.


