# Project Title: Nonstationary ETAS Models for Temporal Variations in Earthquake Occurrences and Swarms in Taiwan

## Overview

This project involves the application of the Epidemic-Type Aftershock Sequence (ETAS) model to analyze seismic data. The project is implemented using of Fortran program, Python, R and the Generic Mapping Tool (GMT) for seismic data analysis.

## Directory Structure


├── data/ │ ├── complete_dataset/ │ │ ├── Taiwan_updated.csv │ │ ├── taiwan_complete.csv │ │ ├── ZMAP_software/ │ │ │ ├── taiwan_mag3.6.csv │ │ │ ├── taiwan_mag1.csv │ ├── Hualien_dataset/ │ │ ├── GDMScatalog.csv │ │ ├── ZMAP_software/ │ │ │ ├── taiwan_data_zmap_2.3.csv │ │ │ ├── taiwan_data_zmap.csv ├── Fortran_prgram/ │ ├── SASeis2006/ │ │ ├── FORTsources/ │ │ │ ├── etas.f │ │ │ ├── etasim.f │ │ │ ├── retas.f │ │ ├── input_annotation/ │ ├── Section3/ │ │ ├── multistage_changepoints/ │ │ │ ├── subperiod_1/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.etas │ │ │ │ ├── output_file/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.res │ │ │ │ ├── etas.error │ │ │ ├── subperiod_2/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.etas │ │ │ │ ├── output_file/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.res │ │ │ │ ├── etas.error │ │ │ ├── subperiod_3/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.etas │ │ │ │ ├── output_file/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.res │ │ │ │ ├── etas.error │ │ ├── single_changepoint/ │ │ │ ├── subperiod_1/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.etas │ │ │ │ ├── output_file/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.res │ │ │ │ ├── etas.error │ │ │ ├── subperiod_2/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.etas │ │ │ │ ├── output_file/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.res │ │ │ │ ├── etas.error │ │ ├── stationary_ETAS/ │ │ │ ├── work.etas │ │ │ ├── etas.open │ │ │ ├── etasim.open │ │ │ ├── output_file/ │ │ │ │ ├── etas.open │ │ │ │ ├── work.res │ │ │ │ ├── simulation_file.zip │ ├── Section4/ │ │ ├── vary_mu/ │ │ │ ├── etas.open │ │ │ ├── work.etas │ │ │ ├── output_file/ │ │ │ │ ├── mu_values_1day.csv │ │ │ │ ├── mu_values_5days.csv │ │ │ │ ├── mu_values_10days.csv │ │ ├── vary_mu_p/ │ │ │ ├── etas.open │ │ │ ├── work.etas │ │ │ ├── output_file/ │ │ │ │ ├── mu_p_values_1day.csv │ │ │ │ ├── mu_p_values_5days.csv │ │ │ │ ├── mu_p_values_10days.csv │ │ │ │ ├── simulation1_same_mag.zip │ │ ├── etas.f │ │ ├── etasim.f │ │ ├── retas.f ├── Python/ │ ├── Section3_changepoints/ │ │ ├── gfortran.ipynb │ │ ├── simulation_lrts.ipynb │ │ ├── decision_tree_output.ipynb │ ├── Section4_moving_windows/ │ │ ├── mu_p_ks1_test.ipynb │ │ ├── J3_statistics_twosample.ipynb ├── R/ │ ├── section2and3plots.R │ ├── log_intensity.R │ ├── monte_carlo_simulation.R │ ├── section4_moving_windows.R ├── GMT/ │ ├── tmp.focal │ ├── tmp.sh │ ├── CN-faults.gmt │ ├── output_png/ │ │ ├── FOCAL.png └── README.md

## Data

- **`complete_dataset/`**:
  - **`Taiwan_updated.csv`**: 
    - The original dataset containing seismic events with magnitude \(M \geq 1\).
  - **`taiwan_complete.csv`**: 
    - A filtered dataset of the original dataset, where events with magnitude \(M \geq 3.6\) are retained.
  - **`ZMAP_software/`**:
    - **`taiwan_mag3.6.csv`** and **`taiwan_mag1.csv`**: 
      - Restructured datasets prepared compatible with ZMAP.
- **`Hualien_dataset/`**:
  - **`GDMScatalog.csv`**: 
    - Dataset of Hualien seismicity.
  - **`ZMAP_software/`**:
    - **`taiwan_data_zmap_2.3.csv`** and **`taiwan_data_zmap.csv`**: 
      - Restructured datasets prepared compatible with ZMAP, threshold magnitude 2.3 and threshold magnitude 1.

## Fortran Program
The core computational tasks of this project are carried out using Fortran Program, specifically using the SASeis2006 package developed by [Ogata (2006)](https://www.ism.ac.jp/~ogata/Ssg/ssg_softwaresE.html). The Fortran code is responsible for performing critical operations such as maximum likelihood estimation (MLE) of ETAS parameters, calculating standard errors, and simulating earthquake sequences using the thinning method.

- **`SASeis2006/`**: 
  - Contains the primary source code for the ETAS model implementation.
  - **`FORTsources/`**: 
    - **`etas.f`**: Computes Maximum Likelihood Estimates (MLE) of ETAS parameters, calculates standard errors, and outputs AIC and log-likelihood values.
    - **`etasim.f`**: Simulates earthquake sequences using the thinning method.
    - **`retas.f`**: Output the residuals based on input file etas.open and work.etas.
  - **`input_annotation/`**: 
    - Annotation files related to the input parameters used in the Fortran programs.
    
- **`Section3/`**:
  - Analysis related to stationary, single change points and multistage ETAS models.
  - **`multistage_changepoints/`** and **`single_changepoint/`**:
    - Each of these directories is divided into different subperiods, with input file etas.open, etasim.open, work.etas, etasim.f and etas.f.
    - **`output_file/`** within each subperiod contains:
      - **`etas.open`**: The output file containing the MLE parameters from the ETAS model.
      - **`work.res`**: Residuals (transformed time) from the model fit.
      - **`etas.error`**: Standard errors with MLE of parameters.

- **`Section4/`**:
  - Analyses focused on the temporal variations of the background rate (\(\mu\)) and aftershock decay rate (\(p\)).
  - **`vary_mu/`** and **`vary_mu_p/`**:
    - The ETAS input files and output data that explore the variations in \(\mu\) and \(p\) over different time intervals.
    - **`output_file/`** within each directory includes CSV files representing the calculated values over different time windows (e.g., 1 day, 5 days, 10 days) and simulation results based on the same magnitude sequence of dataset.

## Python

The Python code in this project are used for executing the Fortran programs, analyzing the output data, detecting change points, and running simulations.

### Directory Breakdown:

- **`Section3_changepoints/`**:
  - Detecting change points in the ETAS model and evaluating their significance.
  - **`gfortran.ipynb`**: 
    - Executes the Fortran programs, processes the outputs, and iteratively computes the combined Akaike Information Criterion (AIC) values at each time point to detect potential change points.
  - **`simulation_lrts.ipynb`**: 
    - Conducts simulations to assess the significance of the detected change points using the Likelihood Ratio Test Statistic (LRTS).
  - **`decision_tree_output.ipynb`**: 
    - Decision tree plot from the output.

- **`Section4_moving_windows/`**:
  - Temporal variations in seismic parameters using a moving window approach.
  - **`mu_p_ks1_test.ipynb`**: 
    - Executes the Fortran programs and use moving windows approach to compare AIC values of different sizes of moving windows. Then run a Kolmogorov-Smirnov (KS) test to test the Poisson hypothesis.
  - **`J3_statistics_twosample.ipynb`**: 
    - Detect changes in the background rate (\(\mu\)) and aftershock decay rate (\(p\)) over time from the outputs of Fortran programs.

## R

R is used for generation of plots that represent the results of the ETAS model (Cumulative number of events versus ordinary and transformed time) and other analyses.
### Directory Breakdown:

- **`section2and3plots.R`**: 
  - This script generates spatial distribution plots of the study region, the cumulative event plots versus ordinary and transformed time. These visualizations are based on the outputs from the Fortran programs.
- **`log_intensity.R`**: 
  - Computes the log-intensity of the 2018 earthquake sequence.
- **`monte_carlo_simulation.R`**: 
  - Implements Monte Carlo simulations to estimate the penalty term \( q(N) \) in the model, which is used for evaluating the additional complexity for searching change points.
- **`section4_moving_windows.R`**: 
  - Plots the results of the temporal variations in \(\mu\) and \(p\) as detected by the moving windows analysis.

## GMT

The GMT (Generic Mapping Tool) for visualizing the spatial distribution of Hualien swarms in 2021.

- **`tmp.focal`**: 
  - A file used in GMT for plotting focal mechanisms of earthquakes.
- **`tmp.sh`**: 
  - A shell script that automates the execution of GMT commands.
- **`CN-faults.gmt`**: 
  - Contains the fault line data from GMT Geospatial Data Collection.
- **`output_png/`**: 
  - **`FOCAL.png`**: Output picture.



## References
[1] Ogata, Y. (2006). SASeis2006: https://www.ism.ac.jp/~ogata/Ssg/ssg_softwaresE.html.

[2] Wiemer, S. (2001). A software package to analyze seismicity: Zmap. Seismological Research Letters, 72(3):373–382.

[3] GMT. China geospatial data collection. https://github.com/gmt-china/
china-geospatial-data. Accessed: 2024-08-12.
