# EEG pipeline for automatic preprocessing and feature extraction

This is a workflow that automatically preprocess and extracts features of EEG data. It is designed for resting state data but can also be applied to Evoked Related Potentials. The input data needs to be raw EEG data in BIDS format. The output data are the preprocessed EEG data, the below mentioned brain features, and a visualization of the results in form of a PDF report.

## Description
The preprocessing of the EEG data is based on Pernet et. al 2019. The main steps are:
- **Detection and rejection of bad channels**
- **Re-referencing to the average reference**
- **Independent Component Analysis:** automatic detection and removal of artifactual independent components
- **Interpolation of rejected bad channels**
- **Detection and rejection of bad time segments**

The following brain features are extracted and plotted:
- **Power spectrum (sensor space)**
- **Alpha Peak Frequency (sensor space)**
- **Average power spectrum (source space)**
- **Functional connectivy (source space):** debiased weighted Phase Lag Index and Amplitude Envelope Correlation 
- **Network characterization (source space):** two local graph theory measures (degree, clustering coefficient) and three global measures (global clustering coefficient, global efficiency, smallworldness)

The brain features in the source space are computed separately for four frequency bands of interest (theta, alpha, beta and gamma). 

## Getting started
### Dependencies
* Matlab with Signal Processing and Statistical and Machine Learning Toolboxes (Developed and tested on R2020b; also tested in R2022a)
* EEGLab (Developed and tested with v2022.0) The following EEGlab plugins need to be installed:
    * [bids-matlb-tool](https://github.com/sccn/bids-matlab-tools) v6.1 
    * [bva-io](https://github.com/arnodelorme/bva-io) if your data is in BVA format v1.7
    * [firfilt](https://widmann/firfilt) v2.4
    * [clean_rawdata](https://github.com/sccn/clean_rawdata) v2.7
    * [ICLabel](https://github.com/sccn/ICLabel) v1.3
    * [dipfilt](https://github.com/sccn/dipfit) v4.3
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
* By default the output of the pipeline will be found in a subfolder of your ras data folder called 'derivatives_vYYYY_MM_DD'

## Dataset to test the pipeline
A small resting-state dataset with raw data in BIDS format to test this pipeline is available in: XXXX
Additionally you can also find there the generated files after running the pipeline.

## Authors
Cristina Gil Ávila, cristina.gil@tum.de

## Help
For bugs or problems contact cristina.gil@tum.de or open an issue [here](https://github.com/crisglav/eeg-pipeline)

## License
CC BY-NC-SA
