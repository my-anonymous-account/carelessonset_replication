library("ggplot2")

# make the plots for the separate repository
NA_type <- "all"

# load visualizing functions
source("application/plotfuns.R")

# load results
load(paste0("application/results/johnson2005", 
            NA_type, "_results.Rdata"))

# load LSP and RE scores
load(paste0("application/results/johnson2005", 
            NA_type, "_results_LSP.Rdata"))
load(paste0("application/results/johnson2005", 
            NA_type, "_results_scores.Rdata"))


dir_save <- "application/plots_individual-series"
levels <- c(0.001, 0.005, 0.01)

for(level in levels)
{
  cp <- cp_LS_RE[[as.character(level)]]
  idx_flagged <- which(!is.na(cp))
  
  for(i in seq_along(idx_flagged))
  {
    idx <- idx_flagged[i] 
    pp <- plot_dimension(RE = scores[idx,], LSP = lsp[idx,], onset = cp[idx])
    pp <- pp + ggtitle(paste0("ID=", idx, ", carelessness onset flagged at item ", cp[idx]))
    nam <- paste0("alpha=", level, "/alpha=", level, "_idx=", idx, ".pdf")
    ggsave(filename = nam, plot = pp, device = "pdf", 
           path = dir_save, width = 5, height = 3.5)
  }
}

