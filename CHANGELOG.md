# Changelog

All notable changes to this project will be documented in this file.


## Version 2.0.0

### Added
- Preprocessing event related data is now possible. Parameters 'PreprocEventData','EventMarker', and 'EventBounds' added.
- Log files are created during preprocessing for each recording. A general log file of the preprocessing is also created to check if all the recordings were correctly preprocessed.
- All parameters are stated and described in the README file.
- ICA block (steps 4, 5 and 6) is performed without parallelization if it fails during parallelization.
- Parameter 'RejectBadTimeSegments' was introduced. Now the rejection of bad time segments is optional. If bad segments are not rejected, the continous recording is saved together with a temporal mask indicating the time samples that were marked as bad segments.
- Parameter 'BrainFeatExt' was introduced. Now the extraction of brain features can be ommited.
- Parameter 'Pad' was introduced to select the zero padding during the estimation of the PSD.

### Changed
- ICA block (steps 4, 5 and 6) is performed in a separate function (preprocessing_ICA.m)
- Selection of 'average' ICA block is performed in a separate function (preprocessing_select_ICA_rep.m)
- Parameters have to be defined in a file called params.json. See examples of this file in params_example.json and params_example_preproc_event.json
- Parameter 'DownsamplingRate' is optional now. If empty, downsampling is skipped.
- Parameter 'FuseChanRej' was eliminated.
- Parameters 'EpochLength' and 'EpochOverlap' are optional now. If 'EpochLength' is not specified data is not epoched.
- Plotting preprocessing results is done in batches of 30 recordings.

### Fixed
- Fix bug finding the last recording that was preprocessed.
- Parameter 'Session' is correctly parsed now.


## Version 1.0.0

Initial version.

