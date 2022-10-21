# EEG pipeline for automatic preprocessing and feature extraction

This is a workflow that automatically preprocess and extracts features of EEG data. It is designed for resting state data but can also be applied to Event Related Potentials. The input data needs to be raw EEG data in BIDS format, please check that you comply with the standard [here](https://bids-standard.github.io/bids-validator/). The output data are the preprocessed EEG data, the below mentioned EEG features, and a visualization of the results in form of a PDF report.

## Description
The preprocessing of the EEG data is based on Pernet et. al 2019. The main steps are:

0. Downsampling (optional)
1. Line noise removal
2. High pass filtering and bad channel removal
3. Re-referencing to the average reference
4. Independent Component Analysis automatic detection and removal of artifactual independent components
5. Interpolation of rejected bad channels
6. Detection and rejection of bad time segments
7. Segmentation into epochs 


The following EEG features are extracted and plotted:
1. Power spectrum (sensor space)
2. Alpha Peak Frequency (sensor space)
3. Power topographies (source space)
4. Functional connectivy (source space): phase-based (dwPLI) and amplitude-based (oAEC) 
5. Brain network characteristics (source space): two local graph theory measures (degree, clustering coefficient) and three global measures (global clustering coefficient, global efficiency, smallworldness)

Brain features in the source space are computed separately for four frequency bands of interest (theta, alpha, beta and gamma). 

## Getting started
### Dependencies
* Matlab with Signal Processing and Statistical and Machine Learning Toolboxes (Developed and tested on R2020b; also tested in R2022a)
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
* Download Matlab and the above listed Matlab toolboxes.
* EEGLab
    * Enter credentials for download in https://sccn.ucsd.edu/eeglab/download.php
    * Download zipped EEGLab
    * Unzip downloaded folder and save on local storage
    * Enter path to software files in `define_params.m`
    * Download the aforementioned EEGlab plugins via the EEGLab GUI or in the provided URL. If you choose the second option, add them to the eeglab/plugins folder. 
      Note that in eeglab/plugins there must not be a same-named subfolder for any plugin (e.g. eeglab/plugins/clean_rawdata-master/clean_rawdata is NOT correct;
      it should be eeglab/plugins/clean_rawdata). By default firfilt, clean_rawdata and ICLabel are in the eeglab/plugins folder when downloading EEGlab.
      If not they need to be added.
* Fieldtrip
    * Enter credentials for download in https://www.fieldtriptoolbox.org/download.php
    * Wait for email that includes link to zipped fieldtrip files and download latest version of fieldtrip
    * Unzip downloaded folder and save on local storage
    * Enter path to software files in `define_params.m`
* Brain connectivity toolbox
    * Download toolbox in https://sites.google.com/site/bctnet/
    * Unzip downloaded folder and save on local storage
    * Enter path to software files in `define_params.m`

### Executing the pipeline
* Read the `define_params.m` file and change the settings according to your needs. The minimum parameters that need to be specified are: 
    * Your study name in `params.study`
    * The path to your BIDS raw data in `params.raw_data_path`
    * If not done previously, the paths of the needed toolboxes
* Run `main_pipeline.m`
* By default the output of the pipeline will be found in a subfolder of your raw data folder called 'derivatives_vYYYY_MM_DD'

## Dataset to test the pipeline
A small resting-state dataset with raw data in BIDS format to test this pipeline is available in: XXXX

## Authors
Cristina Gil Ávila, cristina.gil@tum.de

## Help
For bugs or problems contact cristina.gil@tum.de or open an issue [here](https://github.com/crisglav/eeg-pipeline)

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
