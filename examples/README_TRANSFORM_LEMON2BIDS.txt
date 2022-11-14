
# STEPS TO DOWNLOAD LEMON DATASET AND TRANSFORM IT TO BIDS:
# 1. download the raw EEG files and the participants phenotypic file from http://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON/EEG_Raw_BIDS_ID/ and save them in ~/LEMON-MPI/
# 2. Download the demographics file http://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON/Continous_Peripheral_Physiology_During_MRI_MPILMBB_LEMON/Participants_MPILMBB_LEMON.csv and save it in ~/LEMON-MPI/
# 3. In the Unbuntu terminal, 
cd ~/LEMON-MPI/
find -type d -exec rename 's/RSEEG/eeg/' '{}' \; # Rename the subfolders 'RSEEG' to 'eeg'
# 4. Run the script in matlab lemon8min_fieldtrip.m to extract 8 minutes of resting-state eyes closed data (you need to have fieldtrip)
# 5. Run the script in matlab lemon2bids.m to generate the BIDS sidecar files (you need to have fieldtrip)
# 6. Run the pipeline (adjust previously params.m)

