This repository contains the code for the paper "A novel approach for predicting risk of vector-borne disease establishment in marginal temperate environments under climate change: West Nile virus in the UK".

The process to run the code from start to finish to generate risk maps is as follows:

1. Download the climate projections data from the CEDA data repository as this is owned by the Metoffice (see below)
Citable as:  Met Office Hadley Centre (2018): UKCP18 Regional Projections on a 12km grid over the UK for 1980-2080. Centre for Environmental Data Analysis, 15/01/2021. https://catalogue.ceda.ac.uk/uuid/589211abeb844070a95d061c8cc7f604

2. Run the file "Convert_netcdf_to_csv.R" to convert the chosen netcdf climate projections into a set of .csv files (separate file for each grid square).  Note, for this to work you will need to download minimum and maximum daily temperatures for each desired climate model run for each set of years the WNV model is to be run for.

3. Run the WNV model "WNV_model.f90".  This will require that you download the DDE solver from Thompson and Shampine (https://www.radford.edu/~thompson/ffddes/).  This will create a separate output file for each grid square.

4. Run the R file "Convert_WNVmodel_output_to_dataframes.R" to process the output files to produce a dataframe summarising peak minimum infectious rate (as well as some information about vector-host ratios) for each climate run and each time period.

5. Run the code in "Plot_maps.R" to produce maps as shown in Figures 1, 2, 4, SF1 and SF2.

To reproduce Figure 3 requires that steps 1-3 have been run and then run the file "Produce_fig3.R".

