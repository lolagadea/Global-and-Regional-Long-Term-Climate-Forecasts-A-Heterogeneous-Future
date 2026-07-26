# Global and Regional Long-Term Climate Forecasts: A Heterogeneous Future

## Replication Package

**María Dolores Gadea**  
University of Zaragoza (Spain)

**Jesús Gonzalo**  
Universidad Carlos III de Madrid (Spain)

**International Journal of Forecasting (2026)**

---

## Overview

This repository contains the complete replication package for the paper:

> **Global and Regional Long-Term Climate Forecasts: A Heterogeneous Future**

The package has been reorganized using portable relative paths so that all MATLAB scripts can be executed without modifying local directories.

The repository includes:

- MATLAB source code;
- processed datasets;
- intermediate MATLAB objects required to reproduce the published results;
- documentation describing the execution order;
- references to the original public data sources.

The repository is intended to allow full replication of all empirical results, figures and tables reported in the published paper.

---

# Repository structure

```
01_Data_construction_CRUTEM/
02_Globe/
03_Regional/
Atlas/
Short_run_forecast_real_time/
functions/
functions_for/
Data_sources.txt
README.md
```

## Repository contents

### 01_Data_construction_CRUTEM

Construction of the station database from the CRUTEM temperature records.

### 02_Globe

Global temperature distribution analysis, forecast competition, density forecasts, synthetic-control exercises and long-term projections.

### 03_Regional

Regional temperature distribution analysis and long-term forecasting for all geographical regions considered in the paper.

### Atlas

Comparison between the statistical forecasts developed in this paper and the CMIP/IPCC Atlas climate projections.

### Short_run_forecast_real_time

Real-time short-run forecasting exercises complementing the long-horizon analysis.

### functions

General MATLAB utility functions used throughout the project.

### functions_for

Legacy functions preserved to reproduce exactly the forecasting procedures and model-selection results reported in the paper.

---

# Software requirements

The replication package requires:

- MATLAB
- Statistics and Machine Learning Toolbox
- NetCDF support (for Atlas preprocessing)
- LaTeX (optional, only for compiling generated tables)

---

# Data availability

The repository contains all MATLAB codes, processed datasets and intermediate files required to reproduce the empirical results, figures and tables reported in the published paper.

The original climate databases are **not redistributed** because:

- they are publicly available;
- they are maintained by their original providers;
- their size is very large.

The file

```
Data_sources.txt
```

contains the official download links for all original datasets.

---

# External dependencies

The ARFIMA forecasting exercises rely on third-party MATLAB software, including:

- arfima_est_v2
- MFEToolbox

These packages are **not redistributed** in this repository because they are external software developed and maintained by their original authors.

These external packages affect only the ARFIMA forecasting exercises.
All remaining results reported in the paper can be reproduced directly from the material provided in this repository.

Users wishing to reproduce the ARFIMA forecasting results should obtain these packages from their original sources and place them in the appropriate MATLAB search path before executing the corresponding scripts.

---

# General execution order

The recommended execution sequence is:

1. Construct the CRUTEM station database.
2. Run the Globe preliminary analysis.
3. Run the Globe forecast competition.
4. Generate point and density forecasts.
5. Run the Regional analysis.
6. Run the Atlas comparison.
7. Run the Short-run real-time forecasting exercises.

Intermediate MATLAB objects are provided to facilitate exact replication while substantially reducing computational time.

---

# Reproducibility notes

The repository preserves the computational workflow used in the accepted version of the paper.

Minor modifications have been introduced exclusively to:

- remove hard-coded local paths;
- improve portability across operating systems;
- simplify the replication process.

The scientific results are unchanged.

---

# Citation

If you use this material, please cite:

> Gadea, M.D. and Gonzalo, J. (2026). *Global and Regional Long-Term Climate Forecasts: A Heterogeneous Future*. International Journal of Forecasting.

---

# Contact

For questions regarding the replication package, please contact the corresponding authors.

---

## Notes

This repository contains all materials necessary to reproduce the published results **except**:

- third-party external software (ARFIMA toolboxes);
- publicly available raw climate databases.

Both can be obtained directly from their original providers using the information supplied in **Data_sources.txt**.
