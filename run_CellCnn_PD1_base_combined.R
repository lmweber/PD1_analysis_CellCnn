###################################################################################
# Script to run CellCnn: PD-1 melanoma skin cancer data set (collaboration with   #
# Carsten Krieg)                                                                  #
#                                                                                 #
# data set: "baseline data sets 23 and 29 (combined)"                             #
#                                                                                 #
# note there are two data sets: referred to as "baseline" data sets "23" and "29" #
#                                                                                 #
# Lukas Weber, March 2017                                                         #
###################################################################################


# note: need to run from command line with 'Rscript <filename>.R'


library(flowCore)
library(readxl)
library(dplyr)
library(limma)



# ----------------
# load data into R
# ----------------

# data set: "baseline data sets 23 and 29 (combined)"


# load data from .fcs files
files_base23 <- list.files("../data/PD-1 project/CK_2016-06-23_03all/010_cleanfcs", 
                           pattern = "^BASE_CK_2016-06-23_03_[NR]+[1-5]\\.fcs$", 
                           full.names = TRUE)
files_base29 <- list.files("../data/PD-1 project/CK_2016-06-29_03all3/010_cleanfcs", 
                           pattern = "^BASE_CK_2016-06-29-03all_null_[NR]+[0-9]+\\.fcs$", 
                           full.names = TRUE)

files <- c(files_base23, files_base29)
files

data <- lapply(files, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# check column names (note: exclude columns 58 and 59, which contain "beadDist" and "Time")
ix <- 1:57
check_cols <- lapply(data, function(d) pData(parameters(d))$name)
all(sapply(check_cols, function(ch) all(ch[ix]==check_cols[[1]][ix])))


# load panel details from .xlsx spreadsheet
panel <- read_excel("../data/PD-1 project/CK_panels/panel3_v3.xlsx")
panel



# ------------------------------------
# transform data and export .fcs files
# ------------------------------------

# including CD45


# update 'panel' table to include CD45 in CellCnn analysis
panel$transform[panel$fcs_colname == "Y89Di"] <- 1
panel

# marker columns (to tranform and use for CellCnn analysis)
marker_cols <- as.logical(panel$transform)
marker_cols

panel[marker_cols, ]


# match columns using metal names (since .fcs columns are not in same order as in panels spreadsheet)
marker_metals <- panel[marker_cols, ]$fcs_colname
marker_names <- panel[marker_cols, ]$Antigen

markers_ix <- match(marker_metals, pData(parameters(data[[1]]))$name)

# check
all(panel[marker_cols, ]$fcs_colname == unname(colnames(exprs(data[[1]]))[markers_ix]))


# apply 'asinh' transform with cofactor = 5
cofactor <- 5

data <- lapply(data, function(d) {
  e <- exprs(d)
  e[, markers_ix] <- asinh(e[, markers_ix] / cofactor)
  colnames(e)[markers_ix] <- marker_names
  e
})


# export transformed data
for (i in 1:length(data)) {
  filename <- paste0("../data_transformed/", gsub("\\.fcs$", "", basename(files[i])), "_transf.fcs")
  write.FCS(flowFrame(data[[i]]), filename)
}



# -------------------------------------------------------------------------
# generate .csv files with input arguments for CellCnn (in required format)
# -------------------------------------------------------------------------

files_transf <- list.files("../data_transformed", full.names = TRUE)


# vector of condition IDs
condition <- gsub("[0-9]+_transf\\.fcs$", "", 
                  gsub("^BASE_CK_2016-06-2(3|9)(_|-)(03)(all_null_|_)", "", 
                       basename(files_transf)))
condition

# vector of batch (data set) IDs
batch <- gsub("(_|-)03(all_null_|_)(NR|R)[0-9]+_transf\\.fcs$", "", 
              gsub("^BASE_CK_2016-06-", "", 
                   basename(files_transf)))
batch

# vector of sample IDs
samples <- gsub("_transf\\.fcs$", "", 
                gsub("^BASE_CK_2016-06-2(3|9)(_|-)(03)(all_null_|_)", "", 
                     basename(files_transf)))
samples


# create data frame of sample names and conditions (for CellCnn input .csv file)

label <- as.numeric(as.factor(condition)) - 1
label

df_samples <- data.frame(fcs_filename = basename(files_transf), 
                         label = label)
df_samples


# create data frame of column names (markers) (for CellCnn input .csv file)

df_markers <- t(data.frame(marker_names))
df_markers


# save as .csv files

write.csv(df_samples, "../inputs/input_samples.csv", quote = FALSE, row.names = FALSE)

# need to use 'write.table' to allow removing column names
write.table(df_markers, "../inputs/input_markers.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)



# -------------------------
# investigate batch effects
# -------------------------

# consider each data set ("23" and "29") to be a batch

# note: CellCnn cannot accept covariates to deal with batch effects


# summarize data for MDS plot: median marker expression per sample
n_cells <- sapply(data, nrow)
n_cells

smp <- rep(samples, n_cells)
smp <- factor(smp, levels = samples)

# marker columns only
data_MDS <- lapply(data, function(e) {
  e_markers <- e[, markers_ix]
  e_markers
})

data_MDS <- do.call("rbind", data_MDS)

data_MDS <- data.frame(data_MDS, sample = smp)

data_MDS %>% 
  group_by(sample) %>% 
  summarize_all(median) -> 
  df_meds

# rearrange data frame for MDS plot
df_plot <- t(df_meds[, -1])
colnames(df_plot) <- df_meds$sample


# MDS plot: color by condition

pdf("../plots/MDS_plot_condition.pdf", width = 6.5, height = 6.5)

cols_cnd <- as.character(factor(condition, labels = c("blue", "orange")))
plotMDS(df_plot, top = 2000, col = cols_cnd, 
        main = "MDS plot: color by condition (NR vs. R)")
legend("topleft", legend = c("Non responder (NR)", "Responder (R)"), 
       col = c("blue", "orange"), pch = 16)

dev.off()


# MDS plot: color by batch

pdf("../plots/MDS_plot_batch.pdf", width = 6.5, height = 6.5)

cols_bch <- as.character(factor(batch, labels = c("deepskyblue1", "red")))
plotMDS(df_plot, top = 2000, col = cols_bch, 
        main = "MDS plot: color by batch (dataset 23 vs. dataset 29)")
legend("topleft", legend = c("dataset base_23", "dataset base_29"), 
       col = c("deepskyblue1", "red"), pch = 16)

dev.off()



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


