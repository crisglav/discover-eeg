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

The following brain features are extracted:
Sensor space
- Power spectrum
- Peak frequency in the alpha range
Source space
- Power spectrum
- Functional connectivy: debiased weighted Phase Lag Index and Amplitude Envelope Correlation
- Graph theory measures based on the functional connectivity matrices.

### Dependencies
* Matlab. Developed and tested on ML2020b
* Matlab Signal Processing Toolbox and Statistical Toolbox
* EEGLab. Developed and tested with v2022.0
* Fieldtrip. Developed and tested with vXXXX
* Brain connectivity toolbox 2019_03_03 (https://sites.google.com/site/bctnet/)

The following EEGLab plugins need to be downloaded and saved in eeglab/plugins
* bids-matlb-tool v6.1 (https://github.com/sccn/bids-matlab-tools)
* clean_rawdata v2.7 (https://github.com/sccn/clean_rawdata)
* ICLabel v1.3 (https://github.com/sccn/ICLabel)


## Getting started
In define_params.m 
* Define your study name in params.study
* Define the path of your raw data in params.raw_data_path
* Define the path of your preprocessed data and brain features in params.preprocessed_data_path
* Define the path of the stored toolboxes
Run main_pipeline.m 

By default the output t
### Data
Exemplary resting-state data to run this pipeline is available in: 

## Authors
Cristina Gil Avila, cristina.gil@tum.de

## License
CC BY-NC-SA
