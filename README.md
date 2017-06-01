# PD-1 analysis: CellCnn

Code to run CellCnn on PD-1 melanoma skin cancer data set, for collaboration project with Carsten Krieg and Malgorzata Nowicka, UZH.


## CellCnn

For installation details and other information on CellCnn, see: [https://github.com/eiriniar/CellCnn](https://github.com/eiriniar/CellCnn)

Note that we are using CellCnn version as at 24 April 2017 (not version 0.1, which was released later). This version can be installed by reverting to previous commits (i.e. the commit as at 24 April 2017) after cloning the CellCnn repository.


## Contents

R scripts:

- [run_CellCnn_PD1_data23data29_baseline_panel3v3.R](run_CellCnn_PD1_data23data29_baseline_panel3v3.R): Script to run CellCnn for "panel3_v3.xlsx". Data set: "data 23" and "data 29" (combined), baseline only, Non-Responders vs. Responders.

- [results_CellCnn_PD1_data23data29_baseline_panel3v3.R](results_CellCnn_PD1_data23data29_baseline_panel3v3.R): Script to analyze results for "panel3_v3.xlsx": calculate summary statistics, calculate statistical test, generate heatmap.


In progress: Panel 1

- [run_CellCnn_PD1_data23data29_baseline_panel1.R](run_CellCnn_PD1_data23data29_baseline_panel1.R): Script to run CellCnn for "panel1.xlsx". Data set: "data 23" and "data 29" (combined), baseline only, Non-Responders vs. Responders.


