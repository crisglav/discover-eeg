# DISCOVER-EEG: an EEG pipeline for biomarker discovery 
This is a workflow that automatically preprocess, analyzes and visualizes resting state EEG data in Matlab using EEGLab and FieldTrip toolboxes. It has been tested on the [LEMON dataset](https://www.nature.com/articles/sdata2018308), the [TD-BRAIN dataset](https://www.nature.com/articles/s41597-022-01409-z), and the [Chronic Pain EEG dataset](https://osf.io/m45j2/).

The acompany publication in Scientific Data can be found [here](https://www.nature.com/articles/s41597-023-02525-0).

[![DOI](https://zenodo.org/badge/518454921.svg)](https://zenodo.org/badge/latestdoi/518454921)

## Description
![pipeline](https://user-images.githubusercontent.com/18517243/212702747-f03f71de-aaf1-4ffb-81e0-963b8333e22b.jpg)

### Data
The input data needs to be raw EEG data in BIDS format, please check that you comply with the standard with the [BIDS validator](https://bids-standard.github.io/bids-validator/). The output data are the preprocessed EEG data, the below-mentioned EEG features, and a visualization of the preprocessing steps and EEG features.

### Preprocessing
0. Downsampling (optional)
1. Line noise removal
2. High pass filtering and bad channel removal
3. Re-referencing to the average reference
4. Independent Component Analysis automatic detection and removal of artifactual independent components
5. Interpolation of rejected bad channels
6. Detection and rejection of bad time segments
7. Segmentation into epochs 

### Feature extraction
1. Power spectrum (sensor space)
2. Alpha Peak Frequency (sensor space)
3. Power topographies (source space)
4. Functional connectivity (source space): phase-based (dwPLI) and amplitude-based (AEC) 
5. Brain network characteristics (source space): two local graph theory measures (degree, clustering coefficient) and three global measures (global clustering coefficient, global efficiency, smallworldness).

Brain features in the source space are computed separately for four frequency bands of interest (theta, alpha, beta, and gamma). 

## Getting started
### Dependencies
* Matlab with Signal Processing and Statistical Toolbox, Machine Learning Toolbox, and Parallel Computing Toolbox (optional). The code was developed and tested on Matlab R2020b, and also tested in R2022a.
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
* Update the path of the toolboxes in `params_example.json`

### Executing the pipeline
* Read the `params_example.json` file and change the settings according to your needs. The minimum parameters that need to be specified are: 
    * Your study name
    * The path to your BIDS raw data
* Run `main_pipeline.m`
* By default, the output of the pipeline will be found in a subfolder of your raw data folder called 'derivatives_vYYYY_MM_DD'

## Citation
If you use this code in your project, please cite:

Gil Ávila C, Bott FS, Tiemann L, Hohn VD, May ES, Nickel MM, Zebhauser PT, Gross J, Ploner P. DISCOVER-EEG: an open, fully automated EEG pipeline for biomarker discovery in clinical neuroscience. Sci Data 10, 613 (2023). doi:[10.1038/s41597-023-02525-0](https://doi.org/10.1038/s41597-023-02525-0)

Additionally, to cite the specific version of DISCOVER-EEG used in your analyses you can use the following Zenodo reference: [10.5281/zenodo.8207523](https://zenodo.org/record/8207523)

## Help
For bugs or problems contact cristina.gil@tum.de or open an issue [here](https://github.com/crisglav/eeg-pipeline)

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
