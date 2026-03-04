# =============================================================================
# Script name: FigureS5_code_with_twinspan.R
# Author: Estephe Kana, Eugeni Belda, Edi Prifti
# Date created: 2025-01-01
# Purpose: Generate Figure S5 including TWINSPAN method - generalization
#          performance comparison (AUC, Accuracy, F1) across all methods
#          using replica-split validation (ABORE train / MBERE test).
# Inputs:  analyses/analysis_outputs/replicas_analysis_results/all_analysis_results_split_replicas_with_twinspan.Rda
#          analyses/analysis_outputs/rf_analysis_by_replicas/*.Rda
# Outputs: analyses/Figures/FigureS5/FigureS5_Best_update_with_twinspan.pdf
# =============================================================================

# -----------------------------------------------------------------------------
# Commands to run the script:
# cd /home/estephe/Documents/Papers/2026-aime_key_fish_bruvs/Bruvs-analysis-code/analyses/Figures/FigureS5
# Rscript FigureS5_code_with_twinspan.R
# -----------------------------------------------------------------------------

library(predomics)
library(dplyr)
library(purrr)
library(ggplot2)
library(viridis)
library(ggh4x)
library(ggVennDiagram)
library(reshape2)
library(patchwork)
library(tidyr)

# Detect project root automatically (works for any user/machine after cloning)
if (!require("here", quietly = TRUE)) install.packages("here")
base_dir <- here::here()

load(file = file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis", "all_analysis_results_split_replicas_with_twinspan.Rda"))

# ── 1. Best-model evaluation (Predomics only) ─────────────────────────────────

params <- expand.grid(
  i = names(alldatalist.video_split_replicas),
  j = c("prev_0", "prev_10"),
  k = c("bininter", "terinter"),
  l = c("predout.bin", "predout.maxn"),
  m = c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE"),
  stringsAsFactors = FALSE
)

best_model_avaluation <- function(i, j, k, l, m) {
  resultObj  <- alldatalist.video_split_replicas[[i]][[j]][[k]][[l]][[m]]
  best.model <- resultObj$digest$best$model
  X_train    <- resultObj$Comp_data$X_train
  X_test     <- resultObj$Comp_data$X_test
  y_test     <- resultObj$Comp_data$y_test
  clf        <- resultObj$clf
  
  X_test_new <- as.data.frame(matrix(0, nrow = nrow(X_test), ncol = ncol(X_train),
                                     dimnames = list(rownames(X_test), colnames(X_train))))
  common_features <- intersect(best.model$names_, colnames(X_test))
  X_test_new[, common_features] <- X_test[, common_features]
  
  best.model.test <- evaluateModel(
    mod = best.model, X = t(X_test_new), y = y_test, clf = clf,
    eval.all = TRUE, force.re.evaluation = TRUE, mode = "test"
  )
  
  data.frame(
    metrics      = c("auc", "accuracy", "f1"),
    value        = round(c(best.model.test$auc_, best.model.test$accuracy_, best.model.test$f1_), 3),
    comparison   = m,
    source       = ifelse(l == "predout.bin", "presAbs", "maxN"),
    featNum      = length(best.model$names_),
    prev_rate    = j,
    method       = paste("BestModel", k, sep = "_"),
    replica_comp = i,
    stringsAsFactors = FALSE
  )
}

pred_best_models.metrics.df <- pmap_dfr(params, best_model_avaluation)

# ── 2. Load RF metrics ────────────────────────────────────────────────────────

load(file = file.path(base_dir, "analyses", "Figures", "FigureS5", "rf_analysis_by_replicas", "rf_analysis_results_split_by_replicas_with_twinspan.Rda"))
load(file = file.path(base_dir, "analyses", "Figures", "FigureS5", "rf_analysis_by_replicas", "rf_analysis_results_split_by_replicas_BestModel_Predomics.Rda"))

metrics.df <- bind_rows(rf.metrics.df, pred_best_models.metrics.df, rf.metrics.bestModel.df)

# remove NAs
# metrics.df <- na.omit(metrics.df)

metrics.df$method <- factor(
  metrics.df$method,
  levels = c("BestModel_bininter", "BestModel_terinter",
             "RF_bestModel_bininter", "RF_bestModel_terinter",
             "RF_bininter", "RF_Indval", "RF_terinter", "RF_Twinspan")
)

metrics.df$replica_comp <- ifelse(
  metrics.df$replica_comp == "ABORE_vs_MBERE",
  "train on ABORE and\ntest on MBERE",
  "train on MBERE and\ntest on ABORE"
)

rownames(metrics.df) <- NULL

auc.df         <- metrics.df[metrics.df$metrics == "auc",      ]
acc.df         <- metrics.df[metrics.df$metrics == "accuracy",  ]
f1.df          <- metrics.df[metrics.df$metrics == "f1",        ]
auc_prev_10.df <- auc.df[auc.df$prev_rate == "prev_10", ]

p1 <- ggplot(data = auc_prev_10.df, aes(x = method, y = value, fill = featNum)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Models", y = "AUC", fill = "Number of species") +
  scale_fill_viridis(option = "D") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1),
        legend.title = element_text(angle = 90, margin = ggplot2::margin(r = 15, b = 10)),
        strip.text.y = element_text(angle = 90)) +
  facet_nested(replica_comp ~ source + comparison)

# ── 3. Indicator-species tables ───────────────────────────────────────────────

load(file = file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis", "all_analysis_results_split_replicas_BestModels_Predomics.Rda"))

indics_species.df    <- read.csv(file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis", "indics_species_analysis_split_by_replicas_with_twinspan.csv"))
indics_species_BM.df <- read.csv(file.path(base_dir, "analyses", "Figures", "FigureS5", "results_replicas_analysis", "indics_species_analysis_split_by_replicas_BestModels_Predomics.csv"))

all_indics_species.df <- bind_rows(indics_species.df, indics_species_BM.df)

# Harmonise method labels
label_map <- c(terinter           = "Predomics_ter",
               bininter           = "Predomics_bin",
               BestModel_terinter = "BestModel_Predomics_ter",
               BestModel_bininter = "BestModel_Predomics_bin",
               twinspan           = "Twinspan")

all_indics_species.df <- all_indics_species.df %>%
  mutate(across(where(is.character), ~ recode(., !!!label_map)))

# ── 4. Venn diagrams ──────────────────────────────────────────────────────────

generate_venn_plot <- function(venn_data, title_text, set_size = 3) {
  ggVennDiagram(venn_data, set_size = set_size, label_size = 5) +
    ggtitle(title_text) +
    theme(plot.title      = element_text(face = "bold", hjust = 0.5),
          legend.title    = element_text(color = "black", size = 9),
          legend.position = "bottom") +
    guides(fill = guide_legend(title = "IndicSpecies")) +
    scale_fill_distiller(palette = "RdBu", limits = c(0, 50))
}

# Build a single grouped summary instead of nested loops:
# for each unique combination of keys, collect the vector of unique species.
dims <- all_indics_species.df %>%
  distinct(replica_comp, comparison, prev_rate, data, source)

# venn_list_by_method: nested as [[replica_comp]][[comparison]][[prev_rate]][[data]][[source]]
venn_list_by_method <- pmap(dims, function(replica_comp, comparison, prev_rate, data, source) {
  list(
    keys    = list(replica_comp = replica_comp, comparison = comparison,
                   prev_rate = prev_rate, data = data, source = source),
    species = unique(all_indics_species.df$feature[
      all_indics_species.df$replica_comp == replica_comp &
        all_indics_species.df$comparison   == comparison   &
        all_indics_species.df$prev_rate    == prev_rate    &
        all_indics_species.df$data         == data         &
        all_indics_species.df$source       == source
    ])
  )
}) %>%
  # Re-nest into the original 5-level list structure expected downstream
  { setNames(., paste(dims$replica_comp, dims$comparison, dims$prev_rate,
                      dims$data, dims$source, sep = "|")) } %>%
  { lst <- list()
  iwalk(., function(entry, nm) {
    k <- entry$keys
    lst[[k$replica_comp]][[k$comparison]][[k$prev_rate]][[k$data]][[k$source]] <<- entry$species
  })
  lst }

# venn_list_by_source: same data, different nesting order [[replica_comp]][[comparison]][[prev_rate]][[source]][[data]]
venn_list_by_source <- { lst <- list()
walk(dims %>% pmap(list), function(k) {
  lst[[k$replica_comp]][[k$comparison]][[k$prev_rate]][[k$source]][[k$data]] <<-
    venn_list_by_method[[k$replica_comp]][[k$comparison]][[k$prev_rate]][[k$data]][[k$source]]
})
lst }

# Generate venn_by_method plots: one plot per [[replica_comp]][[comparison]][[prev_rate]][[data]]
dims_method <- all_indics_species.df %>% distinct(replica_comp, comparison, prev_rate, data)

venn_by_method <- pmap(dims_method, function(replica_comp, comparison, prev_rate, data) {
  list(
    keys = list(replica_comp = replica_comp, comparison = comparison,
                prev_rate = prev_rate, data = data),
    plot = generate_venn_plot(
      venn_list_by_method[[replica_comp]][[comparison]][[prev_rate]][[data]],
      title_text = paste(replica_comp, comparison, prev_rate, sep = "/"))
  )
}) %>%
  { lst <- list()
  walk(., function(entry) {
    k <- entry$keys
    lst[[k$replica_comp]][[k$comparison]][[k$prev_rate]][[k$data]] <<- entry$plot
  })
  lst }

# Generate venn_by_source plots: one plot per [[replica_comp]][[comparison]][[prev_rate]][[source]]
dims_source <- all_indics_species.df %>% distinct(replica_comp, comparison, prev_rate, source)

venn_by_source <- pmap(dims_source, function(replica_comp, comparison, prev_rate, source) {
  list(
    keys = list(replica_comp = replica_comp, comparison = comparison,
                prev_rate = prev_rate, source = source),
    plot = generate_venn_plot(
      venn_list_by_source[[replica_comp]][[comparison]][[prev_rate]][[source]],
      title_text = paste(replica_comp, comparison, prev_rate, source, sep = "/"),
      set_size = 2.5)
  )
}) %>%
  { lst <- list()
  walk(., function(entry) {
    k <- entry$keys
    lst[[k$replica_comp]][[k$comparison]][[k$prev_rate]][[k$source]] <<- entry$plot
  })
  lst }

# ── 5. Input data for heatmaps ────────────────────────────────────────────────

get_comp_data <- function(replica_comp, prev_rate, comparison, data) {
  predomics_obj <- alldatalist.video_split_replicas[[replica_comp]][[prev_rate]][["terinter"]]
  if (!is.null(predomics_obj)) {
    slot <- if (data == "maxN") "predout.maxn" else "predout.bin"
    predomics_obj[[slot]][[comparison]][["Comp_data"]]
  } else {
    twinspan_obj <- alldatalist.video_split_replicas[[replica_comp]][[prev_rate]][["Twinspan"]]
    slot <- if (data == "maxN") "twinspan.maxn" else "twinspan.bin"
    twinspan_obj[[slot]][[comparison]][["Comp_data"]]
  }
}

dims_input <- all_indics_species.df %>%
  distinct(replica_comp, prev_rate, comparison, data)

input_data <- { lst <- list()
pwalk(dims_input, function(replica_comp, prev_rate, comparison, data) {
  lst[[replica_comp]][[prev_rate]][[comparison]][[data]] <<-
    get_comp_data(replica_comp, prev_rate, comparison, data)
})
lst }

# ── 6. Union of indicator species across all methods ─────────────────────────

dims_union <- all_indics_species.df %>% distinct(data, prev_rate, comparison, replica_comp)

# For each (data, prev_rate, comparison, replica_comp), take the union of
# all method-specific species vectors that live in venn_list_by_method
all_methods_indics_species <- { lst <- list()
pwalk(dims_union, function(data, prev_rate, comparison, replica_comp) {
  lst[[data]][[prev_rate]][[comparison]][[replica_comp]] <<-
    Reduce(union, venn_list_by_method[[replica_comp]][[comparison]][[prev_rate]][[data]])
})
lst }

# ── 7. Heatmap data preparation ───────────────────────────────────────────────

load(file = file.path(base_dir, "data", "video_data_object.rda"))

sm$sample_info$Zone <- recode(sm$sample_info$hab,
                              BAR = "Barrier", BAY = "Bay", LAG = "Lagoon")

dfaims      <- sm$X
dfaims_meta <- sm$sample_info
dfaims_taxo <- sm$db_long %>%
  filter(G_SP %in% rownames(dfaims)) %>%
  distinct(G_SP, Family, Genus, TrophicG) %>%
  tibble::column_to_rownames("G_SP")

dfaims_sqrt <- apply(dfaims, 2, sqrt)

dfaims_sqrt_clust.sp      <- hclust(dist(dfaims_sqrt,    method = "euclidean"), method = "ward.D")
dfaims_sqrt_clust.samples <- hclust(dist(t(dfaims_sqrt), method = "euclidean"), method = "ward.D")

dfaims_sqrt_melt <- as.data.frame(dfaims_sqrt) %>%
  tibble::rownames_to_column("feature") %>%
  melt(id.vars = "feature", variable.name = "variable") %>%
  mutate(
    feature  = factor(feature,  levels = dfaims_sqrt_clust.sp$labels[dfaims_sqrt_clust.sp$order]),
    variable = factor(variable, levels = dfaims_sqrt_clust.samples$labels[dfaims_sqrt_clust.samples$order])
  ) %>%
  left_join(dfaims_meta %>% tibble::rownames_to_column("variable") %>%
              select(variable, hab, Habitat, Zone),
            by = "variable") %>%
  left_join(dfaims_taxo %>% tibble::rownames_to_column("feature") %>%
              select(feature, Family, Genus, TrophicG),
            by = "feature")

heatmap_theme <- theme_bw() +
  theme(strip.text.y = element_text(angle = 45, vjust = 0.5, hjust = 0.7),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(angle = 90, margin = ggplot2::margin(r = 25, b = 10)))

make_heatmap <- function(df) {
  ggplot(df, aes(x = variable, y = feature, fill = value)) +
    geom_tile(colour = NA) +
    scale_fill_gradient(low = "white", high = viridis::viridis(11)) +
    labs(x = "samples", y = "species", fill = "Species Abundance") +
    facet_nested(. ~ Zone + sampling_set, scales = "free_x", space = "free_x") +
    heatmap_theme
}

# All combinations of (data, prev_rate, comparison) for the outer heatmap loop,
# plus replica_comp for the per-replica annotation step
dims_heat_outer <- all_indics_species.df %>% distinct(data, prev_rate, comparison)
dims_heat_full  <- all_indics_species.df %>% distinct(data, prev_rate, comparison, replica_comp)

# Build per-(data, prev_rate, comparison, replica_comp) annotated data frames,
# then collapse to per-(data, prev_rate, comparison) for plotting.
annotated_tiles <- pmap(dims_heat_full, function(data, prev_rate, comparison, replica_comp) {
  comp_data  <- input_data[[replica_comp]][[prev_rate]][[comparison]][[data]]
  training.X <- comp_data$X_train
  testing.X  <- comp_data$X_test
  
  dfaims_sqrt_melt %>%
    mutate(
      data         = data,
      prev_rate    = prev_rate,
      comparison   = comparison,
      replica_comp = replica_comp,
      sampling_set = case_when(
        variable %in% rownames(training.X) ~ gsub("_.*",  "", replica_comp),
        variable %in% rownames(testing.X)  ~ gsub(".*_",  "", replica_comp),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(sampling_set))
})

# All tiles combined (used to filter to indics species below)
all_tiles_df <- bind_rows(annotated_tiles)

# For each outer combination build the two heatmaps
pheat_map.plots     <- list()
pheat_map_all.plots <- list()

pwalk(dims_heat_outer, function(data, prev_rate, comparison) {
  indics_species <- unique(unlist(
    all_methods_indics_species[[data]][[prev_rate]][[comparison]]
  ))
  
  df_all <- all_tiles_df %>%
    filter(data == !!data, prev_rate == !!prev_rate, comparison == !!comparison)
  
  df_ind <- df_all %>% filter(feature %in% indics_species)
  
  pheat_map.plots[[data]][[prev_rate]][[comparison]]     <<- make_heatmap(df_all)
  pheat_map_all.plots[[data]][[prev_rate]][[comparison]] <<- make_heatmap(df_ind)
})

# ── 8. PDF outputs ────────────────────────────────────────────────────────────

# pdf(file = "/data/projects/aime/analyses/bruvs/FigureS5/FigureS5_Best.pdf", h = 13, w = 11)
# hlay <- "
# A
# A
# B
# B
# B
# "
# p1 + ggtitle("A") +
#   pheat_map_all.plots$maxN$prev_10$INSHORE_OFFSHORE + ggtitle("B") +
#   plot_layout(design = hlay)
# dev.off()

pdf(file = file.path(base_dir, "analyses", "Figures", "FigureS5", "FigureS5_Best_update_with_twinspan.pdf"), h = 22, w = 18)
hlay <- "
AAAA
BBCC
BBCC
DDEE
DDEE
"
p1 + ggtitle("A") +
  pheat_map_all.plots$maxN$prev_10$INSHORE_OFFSHORE + ggtitle("B") +
  pheat_map_all.plots$maxN$prev_10$BAR_BAY          + ggtitle("C") +
  pheat_map_all.plots$maxN$prev_10$BAR_LAG          + ggtitle("D") +
  pheat_map_all.plots$maxN$prev_10$BAY_LAG          + ggtitle("E") +
  plot_layout(design = hlay)
dev.off()

