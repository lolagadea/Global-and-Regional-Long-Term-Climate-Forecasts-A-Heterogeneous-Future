This folder contains the raw CRUTEM5 station files used to construct the STATIONS.mat database.

Raw data source:
Climatic Research Unit (CRU), University of East Anglia
CRUTEM5 station data, version 5.0.2.0

Download page:
https://crudata.uea.ac.uk/cru/data/temperature/

The folder CRUTEM.5.0.2.0.station_files contains the raw NetCDF station files downloaded from CRU. These files are read by Reading_folders_files.m to construct STATIONS.mat.

The script detect_crazy_stations.m identifies stations with invalid coordinates or elevation metadata and creates index_crazy.dat.