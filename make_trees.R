if(interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake",
    slots = c(
      input = "list",
      output = "list",
      params = "list",
      wildcards = "list",
      threads = "numeric",
      log = "list",
      resources = "list",
      config = "list",
      rule = "character",
      bench_iteration = "numeric",
      scriptdir = "character",
      source = "function"
    )
  )
  # if manual execution, the paths and parameters need to be adapted
  results_dir <- "~/Documents/tools/CoreSimul/mash_vs_mosaic/mummer_lastz_2/tree_mummer_results_delta_1000_rrho_0.00e+00_tau_1.00e+00_codon/"
  snakemake <- Snakemake(
    input = list(
      fitted_params = list.files(paste0(results_dir, "fitted_params/")),
      inflexion_file = paste0(results_dir,"inflexion_exists.csv"),
      inflexion_percentage = paste0(results_dir, "inflexion_by_cluster.csv")
    ),
    output = list(),
    params = list(
      taxon_csv = "~/Data/simulated_datasets/tree_genomes_delta_1000_rrho_0.00e+00_tau_1.00e+00_codon/taxon.csv",
      cluster_name = "clade",
      fitted_params_dir = paste0(results_dir, "fitted_params/"),
      results_dir = results_dir,
      genome_wise_inflexion = "no",
      genome_lengths = paste0(results_dir, "lengths_distributions/"),
      tree_annotation = "clade",
      filter_min_genomes = 1
        ),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list(),
    config = list(),
    rule = "",
    bench_iteration = 1,
    scriptdir = "",
    source = function(...) {{ wd <- getwd()
      setwd(snakemake@scriptdir)
      source(...)
      setwd(wd) }}
  )
}

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

read_counts <- function(x) {
  x_file <- paste0(x, "_distribution.csv")
  tmp_df <- read_csv(paste0(length_dir, x_file), col_types = cols(col_character(), col_double()))
  tibble("species" = x, "count" = nrow(tmp_df))
}


# Read data --------------------------------------------------------------------

params_dir <- snakemake@params[["fitted_params_dir"]]
taxon_csv <- snakemake@params[["taxon_csv"]]
cluster_name <- snakemake@params[["cluster_name"]]
taxon_df <- read_csv(taxon_csv)
species_list <- taxon_df %>% pull(.data[[cluster_name]]) %>% unique
results_dir <- snakemake@params[["results_dir"]]
length_dir <- snakemake@params[["genome_lengths"]]
inflexions <- read_csv(snakemake@input[["inflexion_file"]])
# the name of the annotation column for the tree
tree_annotation <- snakemake@params[["tree_annotation"]]

fitted_params_files <- list.files(params_dir)

# TODO to parametrize at some point
min_r_infl <- 1

fitted_params <- tibble(
  "log_tau" = numeric(),
  "log_rho" = numeric(),
  "L0" = numeric(),
  "minimum" = numeric(),
  "minimum_minus3" = numeric(),
  "bac_1" = character(),
  "bac_2" = character()
)

filter_min_genomes <- snakemake@params[["filter_min_genomes"]]
species_list <- taxon_df %>%
  group_by(.data[[cluster_name]]) %>%
  filter(n() >= filter_min_genomes) %>%
  pull(.data[[cluster_name]]) %>%
  unique

for(spec_par in fitted_params_files) {
  level_1 <- str_split_1(spec_par, "_vs_")[1]
  level_2 <- str_sub(str_split_1(spec_par, "_vs_")[2], 1, -19)
  tmp_df <- read_csv(file.path(params_dir, spec_par), col_types = cols(.default = col_double()))
  tmp_df <- bind_cols(tmp_df, "bac_1" = level_1, "bac_2" = level_2)
  colnames(tmp_df) <- colnames(fitted_params)
  # réfléchir à dupliquer bac 1 et 2 pour avoir les paires possibles
  # genre bind rows aussi en inversant les bac
  fitted_params <- bind_rows(fitted_params, tmp_df)
}

fitted_params_intra <- fitted_params %>%
  filter(
    grepl(paste(species_list, collapse = "|"), bac_1) &
    grepl(paste(species_list, collapse = "|"), bac_2)
  )

th_comparison_df <- t(combn(species_list, 2, simplify = TRUE))

for (i in seq_len(nrow(th_comparison_df))) {
  select_row <- fitted_params_intra %>%
    filter(
      grepl(th_comparison_df[i, 1], bac_1) & grepl(th_comparison_df[i, 2], bac_2) |
      grepl(th_comparison_df[i, 1], bac_2) & grepl(th_comparison_df[i, 2], bac_1)
    )
  if(nrow(select_row) == 0)
    print(paste("Missing comparison",
                th_comparison_df[i, 1],
                "vs",
                th_comparison_df[i, 2]
                )
          )
}
#
no_inflexion_comps <- inflexions %>%
  filter(infl_exist == "no") %>%
  select(cluster_1, cluster_2)
# construct tree ---------------------------------------------------------------

pseudo_distance <- matrix(NA, nrow = length(species_list), ncol = length(species_list))
colnames(pseudo_distance) <- rownames(pseudo_distance) <- species_list

for (row_i in seq_len(nrow(pseudo_distance))) {
  for (col_i in seq(row_i, ncol(pseudo_distance))) {
    if (col_i == row_i)
      next
    if (inflexions[bacterias_in_reference(inflexions, c(species_list[row_i], species_list[col_i]), "cluster"), "r_infl"] < min_r_infl )
      next

    logtau <- fitted_params_intra %>%
      filter(
        (bacterias_in_reference(., c(species_list[row_i], species_list[col_i])))
      ) %>%
      pull(log_tau)
    pseudo_distance[row_i, col_i] <- logtau
  }
}
# percentage of no_inflexion out of all necessary taus
missing_taus <- 2 * sum(is.na(as.dist(t(pseudo_distance))))/(length(species_list)*(length(species_list)-1))

# filling empty cells with the mean distance over the tree
mean_pseudo_distance <- mean(pseudo_distance, na.rm = TRUE)
for (row_i in seq_len(nrow(pseudo_distance))) {
  for (col_i in seq(row_i, ncol(pseudo_distance))) {
    if(row_i == col_i)
      next
    if (is.na(pseudo_distance[row_i, col_i])) {
      pseudo_distance[row_i, col_i] <- mean_pseudo_distance
    }
  }
}
# to get the proper pseudo distance
tau_distance <- as.dist(t(pseudo_distance))
stopifnot(sum(is.na(tau_distance)) == 0)

# checking treelikeness
delta_res <- delta.plot(10^tau_distance, plot = FALSE)
mean_delta <- mean(delta_res$delta.bar)
write_csv(tibble("mean_delta" = mean_delta, "missing_taus" = missing_taus), paste0(results_dir, "tree_stats.csv"))

png(paste0(results_dir, "delta_plot_hist.png"), width = 10, height = 8, units = "in", res = 300)
delta.plot(10^tau_distance, which = 1)
dev.off()


tree_upgma <- upgma(10^(tau_distance))

# distances and tau/distance comparison ----------------------------------------

coph_distances <- cophenetic(tree_upgma) %>%
  as.data.frame() %>%
  rownames_to_column("bac_1") %>%
  pivot_longer(!bac_1, names_to = "bac_2", values_to = "distance") %>%
  filter(bac_1 != bac_2)


distance_and_fitted <- fitted_params_intra %>%
  anti_join(inflexions %>%
          filter(r_infl < min_r_infl) %>%
        select(cluster_1, cluster_2), by = join_by(bac_1 == cluster_1, bac_2 == cluster_2)) %>%
  inner_join(coph_distances) %>%
  mutate(tau = 10^log_tau) %>%
  mutate(relative_dif = abs(tau - distance)/(tau+distance))

difi <- ggplot(distance_and_fitted, aes(x = distance, y = tau)) +
  geom_point() +
  geom_function(fun = identity)

ggsave(paste0(results_dir, "fitteddistance_vs_founddistance.png"), difi)

difi_hist <- ggplot(distance_and_fitted, aes(x = relative_dif)) +
  geom_histogram(bins = 20, color = "darkblue", fill = "lightblue")

ggsave(paste0(results_dir, "hist_fitteddistance.png"), difi_hist)


# Add taxon, inflexion percentages and counts/inflexion info ------------------------------------------

family_df <- taxon_df %>%
  dplyr::rename(label = all_of(cluster_name))

if (tree_annotation == snakemake@params[["cluster_name"]]) {
  family_df <- family_df %>%
    mutate({{tree_annotation}} := label)
}
counts_df <- family_df %>%
  select(label, genome) %>%
  group_by(label) %>%
  summarise(count = n())

inflexions_per <- read_csv(snakemake@input[["inflexion_percentage"]]) %>%
  rename(label = all_of(cluster_name))

label_order <- tree_upgma %>%
  as_tibble %>%
  filter(!is.na(label)) %>%
  select(label)

fam <- family_df %>%
  distinct(.data[[tree_annotation]], label) %>%
  inner_join(label_order) %>%
  select(label, .data[[tree_annotation]]) %>%
  column_to_rownames("label")


p <- ggtree(tree_upgma) + geom_tiplab()
p <- revts(p) + scale_x_continuous(labels = abs)
    ## scale_x_continuous(labels=function(x) scales::comma(abs(x))) # <-- what you need is actually a function.
  ##

gh <- gheatmap(p, fam,
               colnames = FALSE,
               legend_title = tree_annotation,
               width = 0.1,
               offset = 1e8
               ) +
  scale_x_ggtree() +
  theme_tree2(legend.position = "bottom",
              legend.box = "vertical", legend.margin = margin())

ggsave(paste0(results_dir, tree_annotation, "_tree_big.svg"), gh, width = 8.5, height = 6, dpi = 300)
