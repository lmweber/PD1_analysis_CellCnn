# PD-1 analysis: CellCnn

Code to run CellCnn analyses in the following paper (collaboration project with Carsten Krieg and Malgorzata Nowicka, University of Zurich):

- Krieg C. et al. *High dimensional single cell analysis predicts response to anti-PD-1 immunotherapy*, under review.


## CellCnn

For installation details and other information on CellCnn, see: [https://github.com/eiriniar/CellCnn](https://github.com/eiriniar/CellCnn)


## Contents

A separate R script is included for each analysis. Shell scripts (extension `.sh`) are used to run multiple R scripts.

The R scripts are organized into one directory per panel (panel 1, panel 2, and panel 3v3). Separate scripts are included for analysis of each batch separately ('data 23' and 'data 29') and for both batches combined ('combined'); the 'combined' analyses were used for the final results in the main text.

Each R script contains code to: pre-process data, generate CellCnn inputs, run CellCnn, analyze results, and save output files. Comments at the beginning of each script record which analysis each script corresponds to.


## Main results

The main results referred to in section 'Identification of a monocyte signature using CellCnn' (main text of manuscript) are generated by the following script:

- [analysis_CellCnn_PD1_panel3v3_base_combined.R](https://github.com/lmweber/PD1_analysis_CellCnn/blob/master/panel3v3/analysis_CellCnn_PD1_panel3v3_base_combined.R)


