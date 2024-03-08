# DISCOVER-EEG: an EEG pipeline for biomarker discovery 
This is a workflow that automatically preprocess, analyzes and visualizes resting state EEG data in Matlab using EEGLab and FieldTrip toolboxes. It has been tested on the [LEMON dataset](https://www.nature.com/articles/sdata2018308), the [TD-BRAIN dataset](https://www.nature.com/articles/s41597-022-01409-z), and the [Chronic Pain EEG dataset](https://osf.io/m45j2/).

The accompanying publication in Scientific Data can be found [here](https://www.nature.com/articles/s41597-023-02525-0).

[![DOI](https://zenodo.org/badge/518454921.svg)](https://zenodo.org/badge/latestdoi/518454921)

**New in v.2.0.0. (see CHANGELOG):**
* Preprocessing event-related data is now possible. 
* Log files are created during preprocessing for each recording.
* All parameters are documented in the README.
* Some preprocessing steps and the extraction of brain features are optional now.

## Description
![pipeline](https://user-images.githubusercontent.com/18517243/212702747-f03f71de-aaf1-4ffb-81e0-963b8333e22b.jpg)

### Data
**Input:** the input data needs to be raw EEG data in BIDS format, please check that you comply with the standard with the [BIDS validator](https://bids-standard.github.io/bids-validator/). 

**Output:** the output data are the preprocessed EEG data, the below-mentioned EEG features, a visualization of the preprocessing steps and EEG features, and individual log files for the preprocessing.

### Preprocessing
0. Downsampling (optional)
1. Line noise removal
2. High pass filtering and bad channel removal
3. Re-referencing to the average reference
4. Independent Component Analysis automatic detection and removal of artifactual independent components
5. Interpolation of rejected bad channels
6. Detection and rejection of bad time segments (optional rejection)
7. Segmentation into epochs (optional)

### Feature extraction
1. Power spectrum (sensor space)
2. Alpha Peak Frequency (sensor space)
3. Power topographies (source space)
4. Functional connectivity (source space): phase-based (dwPLI) and amplitude-based (AEC) 
5. Brain network characteristics (source space): two local graph theory measures (degree, clustering coefficient) and three global measures (global clustering coefficient, global efficiency, smallworldness).

By default, brain features in the source space are computed separately for four frequency bands of interest (theta, alpha, beta, and gamma). 

### Preprocessing event-related data
Preprocessing of event-related data is possible since v.1.1.0. Please note that for this type of data, the extraction of brain features is omitted. During event-related data preprocessing, downsampling is not performed (DownsamplingRate: []), bad time segments are not rejected (RejectBadTimeSegments: "off"), and data is not segmented into epochs (EpochLength: []). In addition, the parameters PreprocEventData, EventMarker and EventBounds have to be specified. Please see `params_example_preproc_event.json` for an example.

## Getting started
### Dependencies
* Matlab with Signal Processing Toolbox, Statistics and Machine Learning Toolbox, and Parallel Computing Toolbox (optional). The code was developed and tested on Matlab R2020b, and also tested in R2022a.
* EEGLab (Developed and tested with v2022.0) The following EEGlab plugins need to be installed:
    * [bids-matlb-tool](https://github.com/sccn/bids-matlab-tools) v6.1 
    * [bva-io](https://github.com/arnodelorme/bva-io) if your data is in BVA format v1.7
    * [clean_rawdata](https://github.com/sccn/clean_rawdata) v2.7
    * [cleanline](https://github.com/sccn/cleanline) v2.0
    * [dipfilt](https://github.com/sccn/dipfit) v4.3
    * [firfilt](https://widmann/firfilt) v2.4
    * [ICLabel](https://github.com/sccn/ICLabel) v1.3
* Fieldtrip (tested with fieldtrip-20220104)
* Brain connectivity toolbox (Developed and tested with version 2019_03_03) 

### Installation
* Download the pipeline code [here](https://github.com/crisglav/eeg-pipeline)
* Download Matlab and the above-listed Matlab toolboxes
* Download [EEGLab](https://sccn.ucsd.edu/eeglab/download.php)
    * Download the aforementioned EEGlab plugins via the EEGLab GUI or in the provided URL. If you choose the second option, add them to the eeglab/plugins folder. 
      Note that in eeglab/plugins there must not be a same-named subfolder for any plugin (e.g. eeglab/plugins/clean_rawdata-master/clean_rawdata is NOT correct;
      it should be eeglab/plugins/clean_rawdata).
    * Change [EEGLab preferences](https://eeglab.org/tutorials/misc/EEGLAB_option_menu.html) in EEGLab GUI. Make sure that the box 'Keep at most one dataset in memory' is selected.
* Download [FieldTrip](https://www.fieldtriptoolbox.org/download.php)
* Download the [BCT toolbox](https://sites.google.com/site/bctnet/)

### Test data
You can download a resting-state EEG test dataset in BIDS format here: https://osf.io/mru42/

### Executing the pipeline
* Copy the file `params_example.json` and rename it to `params.json`.
* Update the `params.json` according to your needs. See [below](#parameters-description) for a description of all parameters.
* Run `main_pipeline.m`
* By default, the output of the pipeline will be found in a subfolder of your raw data folder called 'derivatives_vYYYY_MM_DD'

## Parameters description

|Parameter  			|Required    | Type     |Description
|-------------------------------|------------|----------|--------------
|NCores     			|Optional    |Number    |Number of cores to run in parallel used in CleanLineNoise and Independent Component Analysis.
|EEGLabPath 			|Required    |String    |Absolute path to EEGLab.
|FieldtripPath 			|Required    |String    |Absolute path to FieldTrip.
|BrainConnectivityToolboxPath	|Required    |String	|Absolute path to the Brain Connectivity Toolbox.
|StudyName  			|Optional    |String    |Name of the study.
|RawDataPath			|Required    |String    |Absolute path to the raw data in BIDS format.
|PreprocessedDataPath		|Optional    |String 	|Absolute path to store the preprocessed data. If not specified, it will be derivatives_timestamp.
|Session    			|Optional    |String or Array|Label(s) of the session(s) to process (e.g., ["sess1", "sess2"]). If set to [] all the sessions are analyzed.
|Run        			|Optional    |Number or Array|Label(s) of the run(s) to process (e.g., [1, 2, 3]). If set to [] all the runs are analyzed.
|Task       			|Optional    |String or Array|Label(s) of the task(s) to process (e.g., ["baseline", "task"]). If set to [] all the tasks are analyzed.
|BidsChanloc			|Required    |String    |"on" keeps channel locations from BIDS electrodes.tsv file, "off" reads channel locations from MNI template (standard_1005.elc).
|NoseDir    			|Optional    |String    |Indicates the coordinate system of electrodes.tsv according to EEGLab nomenclature, e.g. RAS is '+Y'. This parameter is required if BidsChanloc is "on".
|RefCoord			|Optional    |Object	|Coordinates X,Y,Z of the reference electrode in the same system as electrodes.tsv. This parameter is required if BidsChanloc is "on" and AddRefChannel is "on".
|DownsamplingRate		|Optional    |Number 	|Rate in Hz to which the data will be downsampled. Set to [] to skip this preprocessing step.
|FlatLineCriterion		|Required    |Number or String	|Maximum tolerated flatline duration of a channel in seconds before being rejected. Parameter derived from clean_artifacts(). Default: 5. To deactivate, set to "off".
|ChannelCriterion		|Required    |Number or String	|Minimum channel correlation. If a channel is correlated less than this value to an estimate based on other channels, it will be considered abnormal in the given time window. Parameter derived from clean_artifacts(). Range: 0 - 1. Default: 0.8. To deactivate, set to "off".
|LineNoiseCriterion		|Required    |Number or String	|If a channel has more line noise relative to its signal than this value, in standard deviations, based on the total channel population, it is considered abnormal. Parameter derived from clean_artifacts(). Default: 4. To deactivate, set to "off".
|HighPass			|Required    |Array or String   |Transition band for the initial high-pass filter in Hz. Parameter derived from clean_artifacts(). Default: [0.25, 0.75]. To deactivate, set to "off".
|AddRefChannel			|Required    |String	|"on" adds the reference channel to the data after rereferencing to average reference, "off" does not add the reference channel back. Default: "off".
|NICARepetitions		|Optional    |Number    |Number of repetitons over the steps 4. ICA, 5. Channel interpolation and 6. Bad segment rejection. Default: 10.
|ICLabel			|Required    |Array     |2x7 array  of floats indicating the probability thresholds to select a component. Maximum and minimum probability thresholds in the range 0-1 have to be specified for each of the seven IC categories: 'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise' and 'Other'. If thresholds are set to null, the category will not be considered. Parameter derived from pop_icflag(). For further documentation consult this function. By default only ICs whose probability of being 'Muscle' or 'Eye' is higher than 0.8 are selected: [[null, null], [0.8, 1], [0.8, 1], [null, null],[null, null],[null, null],[null, null]]
|BurstCriterion			|Required    |Number or String|Standard deviation cutoff for removal of bursts via ASR. Data portions whose variance is larger than this threshold relative to the calibration data will be removed. Parameter derived from clean_artifacts(). Default: 20. To deactivate, set to "off".
|WindowCriterion		|Required    |Number or String|Criterion for removing time windows that were not repaired completely. This may happen if the artifact in a window was composed of too many simultaneous uncorrelated sources (for example, extreme movements such as jumps). This is the maximum fraction of contaminated channels that are tolerated in the final output data for each considered window. Parameter derived from clean_artifacts(). Range: 0 - 1. Default: 0.25. To deactivate, set to "off".
|WindowCriterionTolerances      |Optional    |Array     |Power tolerances outside of which a channel in the final output data is considered "bad" in standard deviations relative to a robust EEG power distribution. Any time window in the final output which has more than the tolerated fraction of channel with power outside this range will be considered incompletely repaired and removed from the output. Parameter derived from clean_artifacts(). Default: "[-Inf, 0.7]".
|RejectBadTimeSegments		|Optional    |String    |"on" rejects portions of data containing bad time segments, concatenates remaining clean data and includes boundary events at appended timepoints, "off" returns the original data and a time mask in EEG.etc.clean\_sample\_mask where clean samples are marked with 1 and bad samples with 0. Default: "on".
|EpochLength			|Optional    |Number	|Length in seconds in which to segment the preprocessed data. If set to [], segmentation into epochs is skipped. This parameter is required for subsequent analyses.
|EpochOverlap			|Optional    |Number    |Amount of epoch overlap in range 0-1. O is no overlap, 1 is complete overlap. Default: 0.
|PreprocEventData |Optional   |Boolean   |If true, event-related data are extracted and preprocessed based on the parameters "EventMarker" and "EventBounds". If set to true, data should not be downsampled, bad time segments will not be rejected and data will not be cut into epochs. Default: false.
|EventMarker      |Optional   |String    |Name of the marker of the event of interest. Default: [].
|EventBounds      |Optional   |Array     |Uper and lower limits in seconds that determine the time window of interest surrounding the "EventMarker". For example, [-2, 1] selects a window -2 seconds before the marker and 1 second after the marker. Default: [].
|BrainFeatExtr			|Optional    |Boolean   |If true, the pipeline performs brain feature extraction, if false, preprocessing only is performed. Brain feature extraction will only be performed if EpochLength is defined. Default: true.
|FreqRes			|Optional    |Number    |Frequency resolution for estimating the power spectrum in Hz. By default is set to 1/EpochLength. If a lower resolution is specified, zero-padding will be applied to the required epoch length.
|Pad 				|Optional    |Number    |Length in seconds to which the data will be padded during frequency analysis. The padding will determine the spectral resolution. Parameter derived from ft_freqanalysis. Default: [].
|FreqBand 			|Optional    |Object    |Frequency bands in which to divide the power spectrum. Each element of the object must include the name of the frequency band of interest and the specified range in Hz, e.g. "beta": [13, 30]. By default the analized frequency bands are theta (4 - 7.9), alpha (8 - 12.9), beta (13 - 30) and gamma (30.1 - 80), with the limits defined by the [COBIDAS MEEG specification](https://www.nature.com/articles/s41593-020-00709-0).
|Taper				|Optional    |String    |Taper applied during frequency analysis. E.g. 'dpss' or 'hanning'. Parameter derived from ft_freqanalysis. Default: 'dpss'.
|Tapsmofrq			|Optional    |Number    |Amount of spectral smoothing through multitapering. Required if Taper is set to 'dpss'. Parameter derived from ft_freqanalysis. Default: 1.
|HeadModelPath			|Optional    |String    |Path to the head model or name of the template if included in FieldTrip. Default: "standard_bem.mat" contained in FieldTrip.
|SurfaceModelPath		|Optional    |String    |Path to the surface model or name of the model if included in FieldTrip. Default: "surface_white_both.mat" contained in FieldTrip.
|AtlasPath			|Optional    |String    |Path to the source atlas. Default: Schaefer atlas with 100 ROIs contained in the parcellations folder ("parcellations/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv").
|FreqResConnectivity		|Optional    |Number    |Frequency resolution for the computation of the debiased weighted phase lag index. For example, if this parameter is set to 0.5, fourier components of the virtual channel data in the alpha band will be computed at 8, 8.5, 9 ... 12.5 Hz and dwpli will be estimated at these frequencies. Default: 0.5.
|ConnMatrixThreshold		|Optional    |Number    |Threshold to binarize connectivity matrices when estimating  graph-theory measures. E.g. if set to 0.2, 20% of the strongest connections will be kept. Range: 0 - 1. Default: 0.2.

## Citation
If you use this code in your project, please cite:

Gil Ávila C, Bott FS, Tiemann L, Hohn VD, May ES, Nickel MM, Zebhauser PT, Gross J, Ploner P. DISCOVER-EEG: an open, fully automated EEG pipeline for biomarker discovery in clinical neuroscience. Sci Data 10, 613 (2023). doi:[10.1038/s41597-023-02525-0](https://doi.org/10.1038/s41597-023-02525-0)

Additionally, to cite the specific version of DISCOVER-EEG used in your analyses you can use the following Zenodo reference: [10.5281/zenodo.8207523](https://zenodo.org/record/8207523)

## Help
For bugs or problems contact cristina.gil@tum.de or open an issue [here](https://github.com/crisglav/eeg-pipeline)

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
