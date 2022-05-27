# An open, fully-automated EEG pipeline

This is a workflow that preprocess EEG data and extracts several brain features. The input data needs to be in BIDS format.


## Description
The preprocessing is based on Pernet et. al 2019. The main steps are:
- Detection of bad channels
- Re-referencing to the average reference
- ICA. Detection and removal of artifactual independent components
- Detection of bad time segments

The following brain features are extracted
- Power spectrum (sensor space)
- Peak frequency (sensor space)
- Functional connectivy: debiased weighted PLI and Amplitude Envelope Correlation (source space)
- Network measures based on the FC measures.

## Getting started

### Dependencies
* Developed and tested on ML2020b. Statistical Signal Processing Toolbox needed.
* EEGLab v2022.0
* Fieldtrip
* Brain connectivity toolbox v 2019_03_03(https://sites.google.com/site/bctnet/)

The following EEGLab plugins need to be downloaded and saved in eeglab/plugins
* bids-matlb-tool v6.1 (https://github.com/sccn/bids-matlab-tools)
* clean_rawdata v2.7 (https://github.com/sccn/clean_rawdata)
* ICLabel v1.3 (https://github.com/sccn/ICLabel)

### Data dependencies
* Parcellations from Schaeffer atlas (add publication)

### Installing
* How/where to download pipeline
* Download all the above listed EEGLab plugins and copy them into eeglab/plugins folder
Modifications that need to be made to files/folders. Modify the define_params.m file before running. Change the paths to your paths.

### Executing pipeline
* Type on the matlab console ´main_pipeline´ after you have modified the define_params.m

## Authors
Cristina Gil Avila, cristina.gil@tum.de

## License
Add license
