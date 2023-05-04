# SGRN
SGRN: A Cas12a-driven Synthetic Gene Regulatory Network System (code to accompany manuscript)

This repository contains code (and data) to reproduce the results in our manuscript entitled "SGRN: A Cas12a-driven Synthetic Gene Regulatory Network System".

There are two directories - "flow" and "microscopy".

The "flow" directory includes multiple subdirectories - "set_processing", "results_combining", and "results_tables". The "set_processing" directory includes code for processing raw flow cytometry fcs files (which can be obtained here: [WEBSITE]), generating plots and create results tables for each replicate batch of each experiment reported in the manuscript. The "results_combining" directory includes code for merging each replicate set of an experiment into a single unified table of results for that experiment. The "results_tables" directory includes copies of all the unified tables we generated for the study. Finally, the "flow" directory also includes R code ("linear_model_code.R") for reproducing all the results figures and running the linear modeling-based tests of significance.

The "microscopy" directory contains one subdirectory ("data"), which includes the table of imaging measurements used to determine differential nuclear localization of the GFP-fused constructs reported in the manuscript. This directory also includes the R code ("single_cell_imaging_data.R") necessary to generate a violing plot and the various tests of significance reported. The original images used in this study are also available here: [WEBSITE].
