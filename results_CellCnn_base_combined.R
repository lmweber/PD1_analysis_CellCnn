###########################################################################
# Script to calculate summary statistics for results from previous script #
# ("run_CellCnn_PD1_base_combined.R")                                     #
#                                                                         #
# Lukas Weber, March 2017                                                 #
###########################################################################


# note: run after script 'run_CellCnn_PD1_base_combined.R', to get 'samples' and
# 'condition' vectors


library(dplyr)


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

freq <- sapply(res_sel, function(d) sum(d[, 2]))
freq


# create data frame of results
res <- data.frame(
  sample = samples, 
  condition = factor(condition), 
  n_cells = n_cells, 
  n_selected = freq, 
  freq_pct = round(freq / n_cells * 100, 2)
)

res


# results per group
sem <- function(x) sd(x) / sqrt(length(x))

res %>% 
  group_by(condition) %>% 
  summarize(mean(freq_pct), 
            sd(freq_pct), 
            sem(freq_pct), 
            median(freq_pct)) -> 
  res_cond

res_cond


# save results
sink("../results/summary_statistics_selected_filter_by_sample.txt")
res
sink()

sink("../results/summary_statistics_selected_filter_by_group.txt")
res_cond
sink()


