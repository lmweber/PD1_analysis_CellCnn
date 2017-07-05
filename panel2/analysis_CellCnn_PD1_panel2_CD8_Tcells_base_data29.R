##########################################################################################
# Script to run CellCnn analysis: run CellCnn, calculate summary statistics, calculate
# statistical test, generate heatmap, save output files
# 
# Anti-PD-1 melanoma skin cancer data set (collaboration with Carsten Krieg and Malgorzata
# Nowicka, UZH)
# 
# - panel: "panel2.xlsx"
# - data: batches "data 23" and "data 29", baseline, Non-Responders vs. Responders
# 
# Note: run from command line with 'Rscript <filename>.R'
# 
# Lukas Weber, June 2017
##########################################################################################


# this script: batch "data 29" only


library(flowCore)
library(readxl)
library(dplyr)
library(limma)
library(lme4)
library(multcomp)
library(pheatmap)
library(RColorBrewer)



########
# inputs
########

dataset <- "panel2_CD8_Tcells_base_data29"

fn_metadata_23 <- "../../data/PD-1 project/CK_metadata/metadata_23_02.xlsx"
fn_metadata_29 <- "../../data/PD-1 project/CK_metadata/metadata_29_02.xlsx"

path_23 <- "../../data/PD-1 project/CK_2016-06-23_02_CD8_merging2/010_cleanfcs"
path_29 <- "../../data/PD-1 project/CK_2016-06-29_02_CD8_merging/010_cleanfcs"

fn_panel <- "../../data/PD-1 project/CK_panels/panel2CD8.xlsx"

batch_name <- "29"



###############
# load metadata
###############

# load metadata spreadsheets for each data set ("data 23" and "data 29")
metadata_23 <- read_excel(fn_metadata_23)
metadata_29 <- read_excel(fn_metadata_29)

#View(metadata_23)
#View(metadata_29)

ix_keep <- 6:15


# paths
paths <- c(rep(path_23, length(ix_keep)), rep(path_29, length(ix_keep)))


# filenames
files <- c(metadata_23$filename[ix_keep], metadata_29$filename[ix_keep])


# vector of condition IDs
condition <- gsub("^base_", "", c(metadata_23$condition[ix_keep], metadata_29$condition[ix_keep]))
condition

# vector of sample IDs
samples <- gsub("^base_", "", c(metadata_23$shortname[ix_keep], metadata_29$shortname[ix_keep]))
samples

# vector of batch (data set) IDs
batch <- c(rep("23", length(ix_keep)), rep("29", length(ix_keep)))
batch


# check
data.frame(paths, files, condition, samples, batch)



##################
# load data into R
##################

# load data from .fcs files
fn <- paste(paths, files, sep = "/")
fn

data <- lapply(fn, read.FCS, transformation = FALSE, truncate_max_range = FALSE)


# align column names
# - remove "SampleID" and "beadDist" (columns 45 and 59) from data29
ix_remove_samples <- c(rep(FALSE, sum(batch == "23")), rep(TRUE, sum(batch == "29")))
ix_remove_cols <- rep(list(c(45, 59)), sum(batch == "29"))
data[ix_remove_samples] <- mapply(function(d, ix) {
  d[, -ix]
}, data[ix_remove_samples], ix_remove_cols)


# check column names
check_cols <- lapply(data, function(d) pData(parameters(d))$name)
all(sapply(check_cols, function(ch) all(ch == check_cols[[1]])))


# load panel details from .xlsx spreadsheet
panel <- read_excel(fn_panel)
panel



################
# transform data
################

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



###########################
# investigate batch effects
###########################

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


# MDS plot: color by condition and batch

cnd_bch <- paste(condition, batch, sep = "_")
cnd_bch <- factor(cnd_bch, levels = unique(cnd_bch))
cnd_bch

pdf(paste0("../../plots/", dataset, "/MDS_plot_condition_batch.pdf"), width = 7.5, height = 7.5)

pal <- c("deepskyblue1", "blue", "orange", "red")
cols_cnd_bch <- as.character(factor(cnd_bch, labels = pal))
plotMDS(df_plot, top = 2000, col = cols_cnd_bch, 
        main = "MDS plot: \ncondition (NR vs. R) and \nbatch (data base_23 vs. data base_29)")
legend("bottomright", pch = 16, 
       legend = c("Non responder (NR), data base_23", "Responder (R), data base_23", 
                  "Non responder (NR), data base_29", "Responder (R), data base_29"), 
       col = pal)

dev.off()



#################################################
# subset for batch; export transformed data files
#################################################

is_batch <- batch == batch_name

files <- files[is_batch]
data <- data[is_batch]
condition <- condition[is_batch]
samples <- samples[is_batch]


for (i in 1:length(data)) {
  filename <- paste0("../../data_transformed/", dataset, "/", 
                     gsub("\\.fcs$", "_transf.fcs", files[i]))
  write.FCS(flowFrame(data[[i]]), filename)
}



###########################################################################
# generate .csv files with input arguments for CellCnn (in required format)
###########################################################################

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

write.csv(df_samples, paste0("../../inputs/", dataset, "/input_samples.csv"), 
          quote = FALSE, row.names = FALSE)

# need to use 'write.table' to allow removing column names
write.table(df_markers, paste0("../../inputs/", dataset, "/input_markers.csv"), 
            sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)



###############################
# run CellCnn from command line
###############################

# for installation instructions and examples see: https://github.com/eiriniar/CellCnn


DIR_CellCnn <- "../../../../../CyTOF/differential/CellCnn/CellCnn/"


# run main analysis
cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
             paste0("-f ../../inputs/", dataset, "/input_samples.csv"), 
             paste0("-m ../../inputs/", dataset, "/input_markers.csv"), 
             paste0("-i ../../data_transformed/", dataset, "/"), 
             paste0("-o ../../out_CellCnn/", dataset, "/"), 
             "--export_csv --group_a NR --group_b R")

runtime_main <- system.time(
  system(cmd)
)

runtime_main

sink(paste0("../../runtime/", dataset, "/runtime_main.txt"))
runtime_main
sink()


# export selected cells
cmd <- paste("python", paste0(DIR_CellCnn, "cellCnn/run_analysis.py"), 
             paste0("-f ../../inputs/", dataset, "/input_samples.csv"), 
             paste0("-m ../../inputs/", dataset, "/input_markers.csv"), 
             paste0("-i ../../data_transformed/", dataset, "/"), 
             paste0("-o ../../out_CellCnn/", dataset, "/"), 
             "--plot", 
             "--group_a NR --group_b R", 
             "--filter_response_thres 0.3 --load_results --export_selected_cells")

runtime_select <- system.time(
  system(cmd)
)

runtime_select

sink(paste0("../../runtime/", dataset, "/runtime_select.txt"))
runtime_select
sink()



#########################################################
# analysis: load results for selected filters/populations
#########################################################

files_sel <- gsub("_transf\\.fcs$", "_transf_selected_cells.csv", files_transf)

fn_sel <- paste0("../../out_CellCnn/", dataset, "/selected_cells/", files_sel)

# check files are in same order as previously
data.frame(fn_sel, samples)

res_sel <- lapply(fn_sel, read.csv)



#######################################################################
# analysis: calculate summary statistics for selected filter/population
#######################################################################

# note: top filter only (if multiple)


n_cells <- sapply(res_sel, nrow)
n_cells

n_selected <- sapply(res_sel, function(d) sum(d[, 2]))  # top filter only (if multiple)
n_selected


# create data frame of results
res <- data.frame(
  sample = factor(samples), 
  condition = factor(condition), 
  n_cells = n_cells, 
  n_selected = n_selected, 
  proportion = n_selected / n_cells, 
  prop_pct = round(n_selected / n_cells * 100, 2)
)

res


# results per group
sem <- function(x) sd(x) / sqrt(length(x))

res %>% 
  group_by(condition) %>% 
  summarize(mean(prop_pct), 
            sd(prop_pct), 
            sem(prop_pct), 
            median(prop_pct), 
            mad(prop_pct)) -> 
  res_cond

res_cond


# save results
sink(paste0("../../results/", dataset, "/summary_statistics_selected_filter_by_sample.txt"))
res
sink()

sink(paste0("../../results/", dataset, "/summary_statistics_selected_filter_by_group.txt"))
t(res_cond)
sink()



##########################################################################
# analysis: statistical tests using generalized linear mixed models (GLMM)
##########################################################################

# note: top filter only (if multiple)


# define model formula
# (using observation-level random effects model to account for overdispersion, i.e. random
# intercept term for each sample; see Harrison (2015): https://peerj.com/articles/1114/)

formula_glmer <- proportion ~ condition + (1 | sample)


# fit model
fit <- glmer(formula_glmer, weights = n_cells, family = binomial, data = res)

fit
summary(fit)


# hypothesis test and p-value
contrast <- as.matrix(t(data.frame(conditionR = c(0, 1))))
contrast

hyp <- glht(fit, linfct = contrast)

hyp
summary(hyp)

p_val <- summary(hyp)$test$pvalues
p_val


# save results
sink(paste0("../../results/", dataset, "/mixed_model_test_results.txt"))
summary(hyp)
sink()



#####################################################################################
# analysis: heatmap of marker expression for selected filter/population vs. all cells
#####################################################################################

# note: top filter only (if multiple)


# indices of cells in selected filter/population
ix_sel <- lapply(res_sel, function(r) {
  which(r[, 2] == 1)
})

ix_sel


# load transformed data files and subset cells
fn_transf <- paste0("../../data_transformed/", dataset, "/", files_transf)

data <- lapply(fn_transf, function(f) {
  exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
})

data_sel <- mapply(function(f, ix) {
  d <- exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
  d_sel <- d[ix, , drop = FALSE]
}, fn_transf, ix_sel)


# marker columns only
data_markers <- lapply(data, function(d) {
  d[, markers_ix, drop = FALSE]
})

data_sel_markers <- lapply(data_sel, function(d) {
  d[, markers_ix, drop = FALSE]
})

# check
t(sapply(data_markers, colnames))
t(sapply(data_sel_markers, colnames))


# calculate median marker expression (pooled cells)
meds <- apply(do.call(rbind, data_markers), 2, median)
meds_sel <- apply(do.call(rbind, data_sel_markers), 2, median)

meds
meds_sel


# # alternatively: calculate median marker expression by sample (i.e. samples equally weighted)
# meds <- t(sapply(data_markers, function(d) {
#   apply(d, 2, median)
# }))
# 
# meds_sel <- t(sapply(data_sel_markers, function(d) {
#   apply(d, 2, median)
# }))
# 
# rownames(meds) <- samples
# rownames(meds_sel) <- samples
# 
# meds
# meds_sel
# 
# 
# # calculate overall medians (median of sample medians; i.e. we are weighting samples equally)
# meds_group <- apply(meds, 2, median)
# meds_sel_group <- apply(meds_sel, 2, median)
# 
# meds_group
# meds_sel_group


# meds_plot <- rbind(meds_group, meds_sel_group)
# rownames(meds_plot) <- c("all cells", "selected population")


# create heatmap
meds_plot <- rbind(meds, meds_sel)
rownames(meds_plot) <- c("all cells", "selected")

pheatmap(meds_plot, 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(100), 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         main = "Median (transformed) marker expression: all cells vs. selected population", 
         filename = paste0("../../plots/", dataset, "/heatmap_marker_expression.pdf"), 
         width = 12, height = 2.5)



###################
# copy output files
###################

SAVE_DIR <- paste0("../../files_for_Carsten/", dataset)

system(paste0("cp ../../results/", dataset, "/mixed_model_test_results.txt", " ", SAVE_DIR))
system(paste0("cp ../../results/", dataset, "/summary_statistics_selected_filter_by_group.txt", " ", SAVE_DIR))
system(paste0("cp ../../plots/", dataset, "/heatmap_marker_expression.pdf", " ", SAVE_DIR))
system(paste0("cp ../../plots/", dataset, "/MDS_plot_condition_batch.pdf", " ", SAVE_DIR))
system(paste0("cp ../../out_CellCnn/", dataset, "/plots/selected_population_boxplot_filter*", " ", SAVE_DIR))
system(paste0("cp ../../out_CellCnn/", dataset, "/plots/selected_population_distribution_filter*", " ", SAVE_DIR))
system(paste0("cp -r ../../out_CellCnn/", dataset, "/selected_cells", " ", SAVE_DIR))

system("rm Rplots.pdf")


