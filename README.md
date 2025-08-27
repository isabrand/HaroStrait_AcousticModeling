Haro Strait Acoustic Modeling

This repo aims to model the acoustic environment for Haro Strait in order to design a localization algorithm for killer whales using a vertical hydrophone array. It explores methods of localization including TDOA and matched field processing beamforming. It uses bellhop as well as classic mathematical processing for signal arrivals at various receivers along the array.

Folder and Code Organization:

TDOA:

at:

beamforming:
* beamforming_bellhop.ipynb
* beamforming_clean.ipynb
* beamforming_funcs.py
* beamforming_functions.ipynb
* beamforming_functions_k_matrices.ipynb
* beamforming_method_comparison.ipynb
* beamforming_method_comparison_k_matrices.ipynb
* bellhop_model.ipynb
* bellhop_testing.ipynb

environmental_modelling:
* bathymetry
  * bathy_line_off_sanjuans_through_mooring.csv
  * bathy_line_off_sanjuans_through_mooring_processed.csv
  * haro_strait.grd
  * sanjuan_cablepath.csv
* haro0_2019.07.24_2019.12.25.nc
* haro_strait_bathy_and_cable.ipynb
* haro_strait_current_data.csv
* haro_strait_filtered_data.csv
* soundspeed_model.ipynb

helpful_papers:

mooring_design:

S01.wav



Please don't hesitate to reach out with questions or comments!
