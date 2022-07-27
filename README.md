# Automatic pipeline for EEG preprocessing and brain feature extraction

This is a workflow that preprocess EEG data and extracts several brain features. 
The input data needs to be raw EEG data in BIDS format.

## Description
The preprocessing of the EEG data is based on Pernet et. al 2019. The main steps are:
- Detection and rejection of bad channels
- Re-referencing to the average reference
- Independent Component Analysis. Automatic detection and removal of artifactual independent components.
- Interpolation of rejected bad channels
- Detection and rejection of bad time segments

The following brain features are extracted and plotted:
In the sensor space:
- Power spectrum
- Peak frequency in the alpha range
In the source space, for each frequency band:
- Power spectrum estimate
- Functional connectivy: debiased weighted Phase Lag Index and Amplitude Envelope Correlation
- Graph theory measures based on the each functional connectivity measure: degree, clustering coefficient, global clustering coefficient, global efficiency and smallworldness

A report of each dataset containing information and visualization of intermediate preprocessing steps and computed brain features can be found in path/to/raw/data/reports folder.

## Getting started
### Dependencies
* Matlab with Signal Processing Toolbox and Statistical Toolbox (Developed and tested on ML2020b; also tested in ML2022a and ML2021b, in which plotting might not work)
* EEGLab (Developed and tested with v2022.0) The following EEGlab plugins need to be installed:
    * bids-matlb-tool v6.1 (https://github.com/sccn/bids-matlab-tools)
    * (optional) bva-io (https://github.com/arnodelorme/bva-io) is needed if your data is in BVA format
    * firfilt (https://widmann/firfilt)
    * clean_rawdata v2.7 (https://github.com/sccn/clean_rawdata)
    * ICLabel v1.3 (https://github.com/sccn/ICLabel)
* Fieldtrip (Developed and tested with revision ee916f5e5)
* Brain connectivity toolbox (Developed and thested with version 2019_03_03) 

### Installation
* Download the pipeline code in XXXXX
* Download Matlab and the above listed Matlab toolboxes.
* EEGLab
    * Enter credentials for download in https://sccn.ucsd.edu/eeglab/download.php
    * Download zipped EEGLab
    * Unzip donloaded folder and save on local storage
    * Enter path to software files in define_params.m file
    * Download the aforementioned EEGlab plugins via the EEGLab GUI or in the provided URL. If you choose the second option, add them to the eeglab/plugins folder. 
      Note that in eeglab/plugins there must not be a same-named subfolder for any plugin (e.g. eeglab/plugins/clean_rawdata-master/clean_rawdata is NOT correct;
      it should be eeglab/plugins/clean_rawdata). By default firfilt, clean_rawdata and ICLabel are in the eeglab/plugins folder when downloading EEGlab.
      If not they need to be added.
* Fieldtrip
    * Enter credentials for download in https://www.fieldtriptoolbox.org/download.php
    * Wait for email that includes link to zipped fieldtrip files and download latest version of fieldtrip
    * Unzip downloaded folder and save on local storage
    * Enter path to software files in define_params.m file
* Brain connectivity toolbox
    * Download toolbox in https://sites.google.com/site/bctnet/
    * Unzip downloaded folder and save on local storage
    * Enter path to software files in define_params.m file

### Executing the pipeline
* Read the define_params.m file and change the settings according to your needs. The minimum parameters that need to be specified are: 
    * Your study name in params.study
    * The path to your BIDS raw data in params.raw_data_path
    * If not done previously, the path of the needed toolboxes
* Run main_pipeline.m 
* By default the output of the pipeline will be found in a 'derivatives_vYYYY_MM_DD' subfolder of your path/to/raw/data

## Dataset to test the pipeline
A small resting-state dataset with raw data in BIDS format to test this pipeline is available in: XXXX
Additionally you can also find here the derivatives of the pipeline.

## Authors
Cristina Gil Avila, cristina.gil@tum.de

## Help
For bugs or problems contact cristina.gil@tum.de or open an issue in XXXX

## License
CC BY-NC-SA
