## this script generates the plots that we use as examples for reconstruction errors and LSP sequences
# for generating the data in the plots we use the 'main' simulation design, that's why this script is here

rm(list = ls()) ; cat("\014")
library("ggplot2")

######### generate data to make the plots from ########
source("R/contamination_illustration.R")
source("R/dgp_johnson2005_cleandata.R")
source("R/benchmarks.R")
source("R/segmentation_performance_measures.R")


# initialize
R              <- 100L
n              <- 500L
scale_size     <- 10L
num_scales     <- 30L
num_likert     <- 5L
batch_size     <- 10L
scale_sizes    <- rep(10L, 30L) # for the option 'factors' in reliability benchmark method
alpha          <- c(0.01, 0.005, 0.001) 
mc_cores <- 40L

# directory to save results in
dir_save <- "illustration/"
prefix <- "exampleplot_"
suffix_type <- "j05sim-main"
suffix_location <- "mixed"


# specify how many samples become careless
prop_crls_arr <- c(0.02, 0.04, seq(from = 0.1, to = 0.6, by = 0.1), 0.8, 1.0)


# number of items
num_items <- scale_size * num_scales

# number alphas
num_alphas <- length(alpha)

## randomly sample the carelessness cutoffs
# we deliberately sample too many (n instead of num_careless), this makes the later code easier
# set.seed(3124851)
set.seed(20260301)

## randomly sample indices of careful respondents (stays fixed across repetitions)
careless_indices_ls <- sample_careless_idx_types(n = n, prop_crls = prop_crls_arr)
careless_start <- carelessness_start_fn(n = n, type = suffix_location)

a <- 4L

prop_crls <- prop_crls_arr[a]
careless_indices_ls_a <- careless_indices_ls[[a]]
idx_pattern <- careless_indices_ls_a$pattern
idx_random <- careless_indices_ls_a$random
idx_careless <- careless_indices_ls_a$careless
idx_careful <- careless_indices_ls_a$careful
group_size_crls <- careless_indices_ls_a$num_careless_per_group


# sample uncontaminated data (ordered as given by johnson2005)
data_ls <- generate_data_johnson05marginals(n)
data <- data_ls$data

# contaminate data
# NOTE: we make the implicit assumption here that careless respondents will
# always choose an answer category (so no response isn't an option). In simulated data, there are no NAs either
data_random_tmp <- data_contamination(data_uncontam  = data[idx_random,,drop = FALSE],
                                      num_likert     = num_likert, 
                                      careless_start = careless_start[idx_random], 
                                      careless_end   = rep(num_items, group_size_crls), 
                                      type           = "random")
data_pattern_tmp <- data_contamination(data_uncontam  = data[idx_pattern,,drop = FALSE],
                                       num_likert     = num_likert, 
                                       careless_start = careless_start[idx_pattern], 
                                       careless_end   = rep(num_items, group_size_crls), 
                                       type           = "pattern")


# merge with uncontaminated data
data_contam                <- data
data_contam[idx_pattern,]  <- data_pattern_tmp
data_contam[idx_random,]   <- data_random_tmp

## run method on 2-dim series
crssonset <- carelessonset::carelessonset(responses = data_contam, 
                                          num_scales = num_scales,
                                          num_likert = num_likert,
                                          longstring = TRUE,
                                          alpha = alpha,
                                          mc_cores = mc_cores,
                                          epochs = 100, 
                                          seed_tf = 12345, 
                                          batch_size = batch_size)


## save the object thus far
save.image("illustration/example_plots.RData")
# load("illustration/example_plots.RData")

# true onset
onset_true <- careless_start
onset_true[idx_careful] <- NA_integer_ # set onset location of truly careful respondents to NA (expected by fun below)

flagged <- crssonset$changepoints$`0.001`
teststat_cutoff <- 192.5 # cutoff value in Shao & Zhang (2010) for d = 2 and alpha = 0.1%

flaggedpositives <- which(!is.na(flagged) & !is.na(onset_true))
truespositives <- sort(idx_careless)



### inconsistent respondent:

# for randoms
idx_rnd <- 454
flagged[idx_rnd] # 181
onset_true[idx_rnd] # 183

## only plot RE series (w/o onset line)
df_re <- data.frame(RE = crssonset$series$RE[idx_rnd,], x = seq_len(num_items), type = "RE")
p_re <- ggplot(df_re, mapping = aes(x = x, y = RE)) +
  theme_bw() + 
  geom_line() +
  # facet_grid(rows = vars(type)) +
  xlab("Item index") +
  # ylab("Score")
  ylab("RE")

ggsave(filename = paste0(prefix, "inconsistentrespondent_REseries.pdf"), plot = p_re, device = "pdf",
       path = dir_save, width = 0.8 * 7.5, height = 0.8 * 2.75)


## now facet with test statistics and LSP
df_rnd <- rbind(
  data.frame(value = crssonset$series$RE[idx_rnd,], x = seq_len(num_items), statistic = "RE"),
  data.frame(value = crssonset$series$LSP[idx_rnd,], x = seq_len(num_items), statistic = "LSP"),
  data.frame(value = c(0, crssonset$teststatistics$Tn[idx_rnd,]), x = seq_len(num_items), statistic = "Test statistic")
)
df_rnd$type <- "Inconsistent carelessness"
levels_statistics <- c("RE", "LSP", "Test statistic")
df_rnd$statistic <- factor(df_rnd$statistic, levels = levels_statistics)

df_hline <- data.frame(statistic = factor(levels_statistics, levels = levels_statistics),
                       cutoff = c(NA, NA, teststat_cutoff)) #teststat_cutoffs

p_teststat_rnd <- 
  ggplot(df_rnd, mapping = aes(x = x, y = value)) +
  theme_bw() +
  geom_line() +
  facet_grid(rows = vars(statistic), #cols = vars(type), 
             scales = "free_y", switch = "y") +
  geom_vline(xintercept = flagged[idx_rnd], color = "blue") +
  geom_hline(aes(yintercept = cutoff), data = df_hline, color = "darkgray", linetype = "dashed") +
  xlab("Item index") +
  # ylab("Score")
  ylab(NULL) + 
  theme(strip.background.y = element_blank(), 
        strip.placement = "outside", 
        strip.text.y = element_text(size = 11), 
        strip.switch.pad.grid = unit(0, "pt"))

ggsave(filename = paste0(prefix, "inconsistentrespondent.pdf"), plot = p_teststat_rnd, device = "pdf",
       path = dir_save, width = 0.49 * 7.5, height = 0.49 * 7.5)



### invariable respondent:

# for patterns
idx_invar <- 57
flagged[idx_invar] # 245
onset_true[idx_invar] # 245
data_contam[idx_invar,260:300] # 1 3 1 3 ...

## only plot LSP series (w/o onset line)
df_lsp <- data.frame(LSP = crssonset$series$LSP[idx_invar,], x = seq_len(num_items), type = "LSP")
p_lsp <- ggplot(df_lsp, mapping = aes(x = x, y = LSP)) +
  theme_bw() + 
  geom_line() +
  # facet_grid(rows = vars(type)) +
  xlab("Item index") +
  # ylab("Score")
  ylab("LSP")

ggsave(filename = paste0(prefix, "invariablerespondent_LSPseries.pdf"), plot = p_lsp, device = "pdf", 
       path = dir_save, width = 0.8 * 7.5, height = 0.8 * 2.75)


## now facet with test statistics and RE
df_str <- rbind(
  data.frame(value = crssonset$series$RE[idx_invar,], x = seq_len(num_items), statistic = "RE"),
  data.frame(value = crssonset$series$LSP[idx_invar,], x = seq_len(num_items), statistic = "LSP"),
  data.frame(value = c(0, crssonset$teststatistics$Tn[idx_invar,]), x = seq_len(num_items), statistic = "Test statistic")
)
df_str$type <- "Invariable partial carelessness"
levels_statistics <- c("RE", "LSP", "Test statistic")
df_str$statistic <- factor(df_str$statistic, levels = levels_statistics)

df_hline <- data.frame(statistic = factor(levels_statistics, levels = levels_statistics),
                       cutoff = c(NA, NA, teststat_cutoff)) #teststat_cutoffs

p_teststat_str <- 
  ggplot(df_str, mapping = aes(x = x, y = value)) +
  theme_bw() +
  geom_line() +
  facet_grid(rows = vars(statistic), #cols = vars(type), 
             scales = "free_y", switch = "y") +
  geom_vline(xintercept = flagged[idx_invar], color = "blue") +
  geom_hline(aes(yintercept = cutoff), data = df_hline, color = "darkgray", linetype = "dashed") +
  xlab("Item index") +
  # ylab("Score")
  ylab(NULL) + 
  theme(strip.background.y = element_blank(), 
        strip.placement = "outside", 
        strip.text.y = element_text(size = 11), 
        strip.switch.pad.grid = unit(0, "pt"))

ggsave(filename = paste0(prefix, "invariablerespondent.pdf"), plot = p_teststat_str, device = "pdf", 
       path = dir_save, width = 0.49 * 7.5, height = 0.49 * 7.5)

