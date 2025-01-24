library(ggtree)
library(dplyr)
library(tibble)
library(readr)
library(phangorn)
library(tidyr)
library(ggplot2)
library(janitor)
library(purrr)

source("./utils.R")

# load data
args <- commandArgs(trailingOnly = TRUE)

res_df <- read_csv(args[1])
taxon_df <- read_csv(args[2])
aligner <- args[3]
delta <- as.numeric(args[4])
results_dir <- args[5]

logtau_df <- res_df %>%
  ## filter(aligner == aligner) %>%
  ## filter(delta == delta) %>%
  select(log10tau, species_1, species_2)

species <- unique(c(logtau_df$species_1, logtau_df$species_2))

logtau_mat <- matrix(NA, nrow = length(species), ncol = length(species), dimnames = list(species, species))

for (i in seq_len(nrow(logtau_df))) {
  logtau_mat[logtau_df$species_1[i], logtau_df$species_2[i]] <- logtau_df$log10tau[i]
  logtau_mat[logtau_df$species_2[i], logtau_df$species_1[i]] <- logtau_df$log10tau[i]
}

logtau_dist <- as.dist(logtau_mat)

# delta plot
if (nrow(res_df) > 3) {
  delta_res <- delta.plot(10^logtau_dist, plot = FALSE)
  mean_delta <- mean(delta_res$delta.bar)
  write_csv(tibble("mean_delta" = mean_delta, paste0(results_dir, "tree_stats.csv")))

  png(paste0(results_dir, "delta_plot_hist.png"), width = 10, height = 8, units = "in", res = 300)
  delta.plot(10^logtau_dist, which = 1)
  dev.off()
}
# TREE
tree_upgma <- upgma(10^(logtau_dist))
write.tree(tree_upgma, paste0(results_dir, aligner, "_", delta, "_tree.nwk"))

# distance/tree comparison
coph_distances <- cophenetic(tree_upgma) %>%
  as.data.frame() %>%
  rownames_to_column("species_1") %>%
  pivot_longer(!species_1, names_to = "species_2", values_to = "distance") %>%
  filter(species_1 != species_2)


distance_and_fitted <- logtau_df %>%
  mutate(tau = 10^log10tau) %>%
  inner_join(coph_distances, by = c("species_1", "species_2")) %>%
  mutate(relative_dif = abs(tau - distance)/(tau+distance))

difi <- ggplot(distance_and_fitted, aes(x = distance, y = tau)) +
  geom_point() +
  geom_function(fun = identity)

ggsave(paste0(results_dir, "fitteddistance_vs_founddistance.png"), difi)

difi_hist <- ggplot(distance_and_fitted, aes(x = relative_dif)) +
  geom_histogram(bins = 20, color = "darkblue", fill = "lightblue")

ggsave(paste0(results_dir, "hist_fitteddistance.png"), difi_hist)


# tree plot

p <- ggtree(tree_upgma) + geom_tiplab()

ggsave(paste0(results_dir, aligner, "_", delta, "_tree.png"), p, width = 8.5, height = 6, dpi = 300)
