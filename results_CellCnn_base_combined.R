#########################################################################################
# Script to calculate summary statistics, statistical test, heatmap (using results from #
# previous script: "run_CellCnn_PD1_base_combined.R")                                   #
#                                                                                       #
# Lukas April, April 2017                                                               #
#########################################################################################


# note: run previous script 'run_CellCnn_PD1_base_combined.R' to get the following 
# vectors: 'samples', 'condition', 'markers_ix'


library(dplyr)
library(flowCore)
library(lme4)
library(multcomp)
library(pheatmap)
library(RColorBrewer)



# -------------------------------------------
# load results for selected filter/population
# -------------------------------------------

files_sel <- list.files("../out_CellCnn/selected_cells", pattern = "\\.csv$", full.names = TRUE)

# check files are in same order as previously
samples
files_sel

res_sel <- lapply(files_sel, read.csv)



# -----------------------------------------------------------
# calculate summary statistics for selected filter/population
# -----------------------------------------------------------

n_cells <- sapply(res_sel, nrow)
n_cells

n_selected <- sapply(res_sel, function(d) sum(d[, 2]))
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
sink("../results/summary_statistics_selected_filter_by_sample.txt")
res
sink()

sink("../results/summary_statistics_selected_filter_by_group.txt")
t(res_cond)
sink()



# ---------------------------------------------------------------
# statistical tests: using generalized linear mixed models (GLMM)
# ---------------------------------------------------------------

# define model formula
formula_glmer <- proportion ~ condition + (1|sample)


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
sink("../results/mixed_model_test_results.txt")
summary(hyp)
sink()



# -------------------------------------------------------------------------
# heatmap of marker expression for selected filter/population vs. all cells
# -------------------------------------------------------------------------

# indices of cells in selected filter/population
ix_sel <- lapply(res_sel, function(r) {
  which(r[, 2] == 1)
})

ix_sel


# load transformed data files and subset cells
files_transf <- list.files("../data_transformed", full.names = TRUE)

data <- lapply(files_transf, function(f) {
  exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
})

data_sel <- mapply(function(f, ix) {
  d <- exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
  d_sel <- d[ix, , drop = FALSE]
}, files_transf, ix_sel)


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
         filename = "../plots/heatmap_marker_exp.pdf", 
         width = 12, height = 2.5)



