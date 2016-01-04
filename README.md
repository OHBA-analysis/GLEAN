# GLEAN

Group Level Exploratory Analysis of Networks                              
Adam Baker 2015

GLEAN is a MATLAB pipeline for identifying patterns of covariation from M/EEG band-limited power using a Hidden Markov Model (HMM) or Independent Component Analysis (ICA), written for [SPM12](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/). This analysis consists of 3 main stages:

1. Computing band-limited amplitude envelopes for single or multiple frequency bands.
2. Reducing the dimensionality of the data via projection to a low-dimensional subspace (via Principal Component Analysis (PCA) or using a parcellation).
3. Decomposition using the HMM or ICA model, from group-concatenated envelopes.

After the decomposition has been run, post-hoc analyses and results may be computed at the session and group-level.

### Requirements:
- MATLAB: http://www.mathworks.com/products/matlab/
- SPM12: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- FSL:  http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/

To view the contents of this toolbox type:
```matlab
help glean
```
For details on how to set up a new GLEAN analysis type:
```matlab
help glean_setup
```
For details on particular GLEAN settings type:
```matlab
help glean_check
```
To run a GLEAN analysis type:
```matlab
run_glean(GLEAN)
```

#### Example dataset: 
An example dataset may be downloaded from [here](https://www.dropbox.com/s/gcci69z9r4toerd/GLEAN_demo.zip?dl=0)


#### If using this toolbox please consider citing the following:

1. "Fast transient networks in spontaneous human brain activity", Baker el al., eLife, 2014
2. "Spectrally resolved fast transient brain states in electrophysiological data", Vidaurre et al., (in rev) 
3. "A symmetric multivariate leakage correction for MEG connectomes", Colclough et al., NeuroImage, 2015
