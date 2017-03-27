###############################################################################
# Script to run CellCnn on PD-1 melanoma skin cancer data set (Carsten Krieg) #
# Lukas Weber, March 2017                                                     #
###############################################################################


# note: need to run from command line with 'Rscript <filename>.R'


library(flowCore)
library(readxl)



# ----------------
# load data into R
# ----------------

# load data from .fcs files
files <- list.files("../data/PD-1 project/CK_2016-06-23_03all/010_cleanfcs", 
                    pattern = "^BASE_CK_2016-06-23_03_[NR]+[1-5]\\.fcs$", 
                    full.names = TRUE)
files

data <- lapply(files, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# check column names
check_cols <- lapply(data, function(d) pData(parameters(d))$name)
sapply(check_cols, function(ch) all(ch==check_cols[[1]]))


# load panel details from .xlsx spreadsheet
panel <- read_excel("../data/PD-1 project/CK_panels/panel3_v3.xlsx")
panel



# ------------------------------------
# transform data and export .fcs files
# ------------------------------------

# update 'panel' table to include CD45 in CellCnn analysis
panel$transform[panel$fcs_colname == "Y89Di"] <- 1
panel

# marker columns (to tranform and use for CellCnn analysis)
marker_cols <- as.logical(panel$transform)
marker_cols

panel[marker_cols, ]


# match columns using metal names (since .fcs columns are not in same order as in panels spreadsheet)
marker_metals <- panel[marker_cols, ]$fcs_colname

markers_ix <- match(marker_metals, pData(parameters(data[[1]]))$name)

# check
all(panel[marker_cols, ]$fcs_colname == unname(colnames(exprs(data[[1]]))[markers_ix]))


# apply 'asinh' transform with cofactor = 5
cofactor <- 5

data <- lapply(data, function(d) {
  e <- exprs(d)
  e[, markers_ix] <- asinh(e[, markers_ix] / cofactor)
  exprs(d) <- e
})


# export transformed data
for (i in 1:length(data)) {
  filename <- paste0("../data_transformed/", gsub("\\.fcs$", "", basename(files[i])), "_transf.fcs")
  write.FCS(flowFrame(data[[i]]), filename)
}



# -------------------------------------------------------------------------
# generate .csv files with input arguments for CellCnn (in required format)
# -------------------------------------------------------------------------

# data frame for .csv file containing sample names and conditions

files_transf <- list.files("../data_transformed", full.names = TRUE)

condition <- gsub("[1-5]_transf\\.fcs$", "", gsub("^BASE_CK_2016-06-23_03_", "", basename(files_transf)))
condition

label <- as.numeric(as.factor(condition)) - 1
label

df_samples <- data.frame(fcs_filename = basename(files_transf), 
                         label = label)
df_samples


# data frame for .csv file containing column names (metals)

df_metals <- t(data.frame(marker_metals))
df_metals


# save as .csv files

write.csv(df_samples, "../inputs/input_samples.csv", quote = FALSE, row.names = FALSE)

# need to use 'write.table' to allow removing column names
write.table(df_metals, "../inputs/input_markers.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)



# -------------------------------
# run CellCnn (from command line)
# -------------------------------

# for installation instructions and examples see: https://github.com/eiriniar/CellCnn


# run main analysis
cmd <- paste("python ../../CellCnn/cellCnn/run_analysis.py", 
             "-f ../inputs/input_samples.csv", 
             "-m ../inputs/input_markers.csv", 
             "-i ../data_transformed/", 
             "-o ../out_CellCnn", 
             "--max_epochs 15 --nrun 10 --train_perc 0.6 --ncell_pooled 5 10 --plot --export_csv", 
             "--group_a NR --group_b R")

runtime_main <- system.time(
  system(cmd)
)

runtime_main

sink("../runtime/runtime_main.txt")
runtime_main
sink()


# export selected cells
cmd <- paste("python ../../CellCnn/cellCnn/run_analysis.py", 
             "-f ../inputs/input_samples.csv", 
             "-m ../inputs/input_markers.csv", 
             "-i ../data_transformed/", 
             "-o ../out_CellCnn", 
             "--plot", 
             "--group_a NR --group_b R", 
             "--filter_response_thres 0.3 --load_results --export_selected_cells")

runtime_select <- system.time(
  system(cmd)
)

runtime_select

sink("../runtime/runtime_select.txt")
runtime_select
sink()


