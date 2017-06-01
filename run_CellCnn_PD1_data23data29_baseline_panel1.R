##########################################################################################
# Script to run CellCnn
# 
# Anti-PD-1 melanoma skin cancer data set
# (collaboration with Carsten Krieg and Malgorzata Nowicka, UZH)
# 
# - data: "data 23" and "data 29" (combined), baseline only, Non-Responders vs. Responders
# - panel: "panel1.xlsx"
# 
# Lukas Weber, June 2017
##########################################################################################


# note: need to run from command line with 'Rscript <filename>.R'


library(flowCore)
library(readxl)
library(dplyr)
library(limma)



# -------------
# load metadata
# -------------

dataset <- "data23data29_baseline_panel1"


# load metadata spreadsheets for each data set ("data 23" and "data 29")
metadata_23 <- read_excel("../data/PD-1 project/CK_metadata/metadata_23_01.xlsx")
metadata_29 <- read_excel("../data/PD-1 project/CK_metadata/metadata_29_01.xlsx")

#View(metadata_23)
#View(metadata_29)


# paths
paths <- c(rep("../data/PD-1 project/CK_2016-06-23_01/010_cleanfcs", length(6:15)), 
           rep("../data/PD-1 project/CK_2016-06-29_01/010_cleanfcs", length(6:15)))


# filenames
files <- c(metadata_23$filename[6:15], metadata_29$filename[6:15])


# vector of condition IDs
condition <- gsub("^base_", "", c(metadata_23$condition[6:15], metadata_29$condition[6:15]))
condition

# vector of sample IDs
samples <- gsub("^base_", "", c(metadata_23$shortname[6:15], metadata_29$shortname[6:15]))
samples

# vector of batch (data set) IDs
batch <- c(rep("23", length(6:15)), rep("29", length(6:15)))
batch


# check
data.frame(paths, files, condition, samples, batch)



# ----------------
# load data into R
# ----------------

# load data from .fcs files
fn <- paste(paths, files, sep = "/")
fn

data <- lapply(fn, read.FCS, transformation = FALSE, truncate_max_range = FALSE)


# align column names
# - remove "SampleID" and "Time" (columns 45 and 59) from data23
# - remove "beadDist" and "Time" (columns 58 and 59) from data29
ix_remove <- c(rep(list(c(45, 59)), sum(batch == "23")), rep(list(c(58, 59)), sum(batch == "29")))
data <- mapply(function(d, ix) {
  d[, -ix]
}, data, ix_remove)


# check column names
check_cols <- lapply(data, function(d) pData(parameters(d))$name)
all(sapply(check_cols, function(ch) all(ch == check_cols[[1]])))


# load panel details from .xlsx spreadsheet
panel <- read_excel("../data/PD-1 project/CK_panels/panel1.xlsx")
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
  filename <- paste0("../data_transformed/", dataset, "/", 
                     gsub("\\.fcs$", "_transf.fcs", files[i]))
  write.FCS(flowFrame(data[[i]]), filename)
}



# -------------------------------------------------------------------------
# generate .csv files with input arguments for CellCnn (in required format)
# -------------------------------------------------------------------------

files_transf <- gsub("\\.fcs$", "_transf.fcs", files)


# create data frame of sample names and conditions (for CellCnn input .csv file)

label <- as.numeric(as.factor(condition)) - 1
label

df_samples <- data.frame(fcs_filename = files_transf, label = label)
df_samples

# re-arrange alphabetically (otherwise CellCnn reads input files in incorrect order)
df_samples <- df_samples[order(df_samples$fcs_filename), ]
df_samples


# create data frame of column names (markers) (for CellCnn input .csv file)

df_markers <- t(data.frame(marker_names))
df_markers


# save as .csv files

write.csv(df_samples, paste0("../inputs/", dataset, "/input_samples.csv"), 
          quote = FALSE, row.names = FALSE)

# need to use 'write.table' to allow removing column names
write.table(df_markers, paste0("../inputs/", dataset, "/input_markers.csv"), 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)



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

pdf(paste0("../plots/", dataset, "/MDS_plot_condition.pdf"), width = 6.5, height = 6.5)

cols_cnd <- as.character(factor(condition, labels = c("blue", "orange")))
plotMDS(df_plot, top = 2000, col = cols_cnd, 
        main = "MDS plot: color by condition (NR vs. R)")
legend("topright", legend = c("Non responder (NR)", "Responder (R)"), 
       col = c("blue", "orange"), pch = 16)

dev.off()


# MDS plot: color by batch

pdf(paste0("../plots/", dataset, "/MDS_plot_batch.pdf"), width = 6.5, height = 6.5)

cols_bch <- as.character(factor(batch, labels = c("deepskyblue1", "red")))
plotMDS(df_plot, top = 2000, col = cols_bch, 
        main = "MDS plot: color by batch (dataset 23 vs. dataset 29)")
legend("topright", legend = c("dataset base_23", "dataset base_29"), 
       col = c("deepskyblue1", "red"), pch = 16)

dev.off()


# MDS plot: color by condition and batch

cnd_bch <- paste(condition, batch, sep = "_")
cnd_bch <- factor(cnd_bch, levels = unique(cnd_bch))
cnd_bch

pdf(paste0("../plots/", dataset, "/MDS_plot_condition_batch.pdf"), width = 7.5, height = 7.5)

cols_cnd_bch <- as.character(factor(cnd_bch, labels = c("blue", "orange", "deepskyblue1", "red")))
plotMDS(df_plot, top = 2000, col = cols_cnd_bch, 
        main = "MDS plot: color by\ncondition (NR vs. R) and \nbatch (dataset base_23 vs. dataset base_29)")
legend("topright", pch = 16, 
       legend = c("Non responder (NR), dataset base_23", "Responder (R), dataset base_23", 
                  "Non responder (NR), dataset base_29", "Responder (R), dataset base_29"), 
       col = c("blue", "orange", "deepskyblue1", "red"))

dev.off()



# -------------------------------
# run CellCnn (from command line)
# -------------------------------

# for installation instructions and examples see: https://github.com/eiriniar/CellCnn


DIR_CellCnn <- "../../../../CyTOF/differential/CellCnn/CellCnn/"


# run main analysis
cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
             paste0("-f ../inputs/", dataset, "/input_samples.csv"), 
             paste0("-m ../inputs/", dataset, "/input_markers.csv"), 
             paste0("-i ../data_transformed/", dataset, "/"), 
             paste0("-o ../out_CellCnn/", dataset, "/"), 
             "--max_epochs 15 --nrun 10 --train_perc 0.6 --ncell_pooled 5 10 --plot --export_csv", 
             "--group_a NR --group_b R")

runtime_main <- system.time(
  system(cmd)
)

runtime_main

sink(paste0("../runtime/", dataset, "/runtime_main.txt"))
runtime_main
sink()


# export selected cells
cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
             paste0("-f ../inputs/", dataset, "/input_samples.csv"), 
             paste0("-m ../inputs/", dataset, "/input_markers.csv"), 
             paste0("-i ../data_transformed/", dataset, "/"), 
             paste0("-o ../out_CellCnn/", dataset, "/"), 
             "--plot", 
             "--group_a NR --group_b R", 
             "--filter_response_thres 0.3 --load_results --export_selected_cells")

runtime_select <- system.time(
  system(cmd)
)

runtime_select

sink(paste0("../runtime/", dataset, "/runtime_select.txt"))
runtime_select
sink()


