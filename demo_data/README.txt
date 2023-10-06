Overview
-------------
This is a small version of the Chronic pain dataset (https://doi.org/10.17605/OSF.IO/M45J2). 

The dataset contains contains eyes-open and eyes-closed electroencephalographic (EEG) resting-state data of 4 patients with chronic pain and it serves for testing purposes. The complete datasest is part of an ongoing project that plans to study medication effects on chronic pain.
This data was recorded between March 2022 and November 2022 in the Klinikum Rechts der Isar (Munich, Germany) and will be updated in the future including more datasets and the clinical characteristics of the participants.

Citing this dataset
------------------------
Please cite as follows:
Zebhauser, P. T., Gil Ávila, C., & Ploner, M. Chronic pain EEG dataset. https://doi.org/10.17605/OSF.IO/M45J2

License
----------
This dataset is made available under the Creative Commons license type CC BY-NC-SA (https://creativecommons.org/about/cclicenses/). Please note that you may distribute, remix, adapt and build upon the material in any medium or format for noncommercial purposes only and only so long as attribution is given to the creator. If you remix, adapt, or build upon the data, you must license the modified material under identical terms. 
Please email to markus.ploner@tum.de if you plan to use the data for publication.

Format
----------
The dataset is formatted according to the Brain Imaging Data Structure (BIDS) for EEG (Pernet et al., 2019, Sci Data). See the 'dataset_description.json' file for the specific version of BIDS used. Generally, you can find data in .tsv files and descriptions in the accompanying .json files.
The data was recorded using mobile, dry-electrode EEG System CGX-32 (Brain Products, Munich, Germany) and stored in the BrainVision data format consisting of three separate files (a text header file (.vhdr) containing meta data, a text marker file (.vmrk) containing information about events in the data and a binary data file (.eeg) containing the voltage values of the EEG). Details of the BrainVision data format can be found on the Brain Products website (https://www.brainproducts.com/productdetails.php?id=21&tab=5).
This dataset was validated with the BIDS validator v1.9.9