library(reshape2)
library(data.table)
library(predomics)
library(dplyr)
library(ggplot2)

# load Random Forest generalization results
load(file = "~/Downloads/Mysubmissions/aime_07052025/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/all_analysis_results_80_20.Rda")

# rename seed levels
names(alldatalist.video_80_20) <- ifelse(names(alldatalist.video_80_20)=="seed_10", "seed_10", gsub("seed_", "seed_0", names(alldatalist.video_80_20)))

# Define parameters
params <- expand.grid(
  i = names(alldatalist.video_80_20),
  j = c("prev_0", "prev_10"),
  k = c("bininter", "terinter"),
  l = c("predout.bin", "predout.maxn"),
  m = c("BAR_BAY", "BAR_LAG", "BAY_LAG", "INSHORE_OFFSHORE"),
  stringsAsFactors = FALSE
)

# Function to evaluate the best model on unknow data (test data)
best_model_avaluation <- function(row) {
  i <- row$i
  j <- row$j
  k <- row$k
  l <- row$l
  m <- row$m
  
  # Retrieve the predomics results object
  resultObj <- alldatalist.video_80_20[[i]][[j]][[k]][[l]][[m]]
  best.model <- resultObj$digest$best$model
  X_train <- resultObj$Comp_data$X_train
  y_train <- resultObj$Comp_data$y_train
  X_test <- resultObj$Comp_data$X_test
  y_test <- resultObj$Comp_data$y_test
  clf <- resultObj$clf
  
  # Build new X_test with same features as X_train
  X_test_new <- as.data.frame(matrix(0, nrow = nrow(X_test), ncol = ncol(X_train)))
  rownames(X_test_new) <- rownames(X_test)
  colnames(X_test_new) <- colnames(X_train)
  
  common_features <- intersect(best.model$names_, colnames(X_test))
  X_test_new[, common_features] <- X_test[, common_features]
  
  # Test the best model on X_test_new
  best.model.test <- evaluateModel(
    mod = best.model,
    X = t(X_test_new),
    y = y_test,
    clf = clf,
    eval.all = TRUE,
    force.re.evaluation = TRUE,
    mode = "test"
  )
  
  # get metrics (AUC, ACC, f1)
  best.model.auc <- round(best.model.test$auc_, digits = 3)
  best.model.acc <- round(best.model.test$accuracy_, digits = 3)
  best.model.f1 <- round(best.model.test$f1_, digits = 3)
  
  # Store metrics
  data.frame(
    metrics = c("auc", "accuracy", "f1"),
    value = c(best.model.auc, best.model.acc, best.model.f1),
    comparison = m,
    source = ifelse(l == "predout.bin", "presAbs", "maxN"),
    featNum = length(best.model$names_),
    prev_rate = j,
    seed = i,
    method = paste("BestModel", k, sep = "_")
  )
}

# Apply the function to each row of parameter combinations
results <- lapply(seq_len(nrow(params)), function(row_idx) {
  row <- params[row_idx, ]
  best_model_avaluation(row)
})

# bind list into a dafarame
pred_best_models.metrics.df <- bind_rows(results)
pred_best_models.metrics.df$seed= ifelse(pred_best_models.metrics.df$seed == "seed_10", "seed_10", gsub("seed_", "seed_0", pred_best_models.metrics.df$seed))

# load rf aanalysis results on 80_20
load(file = "~/Downloads/Mysubmissions/aime_07052025/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/rf_analysis_results_80_20.Rda")

# bind performances metrics for predomics best models and rf
metrics.list <- list(rf.metrics.df, pred_best_models.metrics.df)
metrics.df <- bind_rows(metrics.list)

# reorder seed from 01, 02, to 10
metrics.df= metrics.df[order(metrics.df$seed), ]

# remove rownames of metrics df
rownames(metrics.df) <- NULL

# filtering
metrics_filtered.df <- metrics.df[metrics.df$metrics=="auc" & metrics.df$prev_rate == "prev_10",]

# split by comparison
# sub_metrics_filtered.list <- split(metrics_filtered.df, metrics_filtered.df$metrics)

# convert the dataframe to wide
compute_mean <- function(subdt, var) {
  mean_subdt.df <- subdt %>%
    group_by(comparison, source, method) %>%
    summarise(mean_value = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  out <- reshape2::dcast(mean_subdt.df, method ~ comparison + source, value.var = "mean_value")
  setnames(out, "method", "group")
  rownames(out) <- out$group  
  out$group <- NULL
  return(out)
}

# compute mean accuracy across seeds
# mean_dcasted.df <- lapply(sub_metrics_filtered.list, compute_mean, var = "value")
mean_dcasted.list <- list()
mean_dcasted.list$AUC <- compute_mean(metrics_filtered.df, var = "value")

# compute mean of feature number by method
mean_dcasted.list$FeatNum <- compute_mean(metrics_filtered.df, var = "featNum")

library(fmsb)
library(gridExtra)

# Define colors
colors_border <- c("#9467bd", "#ff7f0e", "#648C67", "#8195B2", "#CC3837")
colors_in <- scales::alpha(colors_border, 0.15)

# Function to create a radar chart
# make_radar_plot <- function(df, title) {
#   # Add max and min rows required by fmsb
#   df_radar <- rbind(
#     max = rep(max(df)+0.02, ncol(df)),
#     min = rep(0, ncol(df)),
#     df
#   )
#   
#   # Add small random jitter to avoid exact overlap
#   df_radar_jitter <- df_radar
#   df_radar_jitter[-c(1,2), ] <- df_radar_jitter[-c(1,2), ] +
#     matrix(runif(n = (nrow(df_radar)-2) * ncol(df_radar),
#                  min = -0.02, max = 0.02),
#            ncol = ncol(df_radar))
#   
#   # Create radar chart
#   # Capture radar plot as a recordable object
#   p <- recordPlot({
#     radarchart(
#           df_radar_jitter, 
#           pcol = colors_border[1:nrow(df)],   # color per method
#           plwd = 2,
#           # pfcol= colors_in,
#           plty = 1,
#           cglcol = "grey",
#           cglty = 1,
#           axislabcol = "grey",
#           # caxislabels = seq(0,1,0.25),
#           # caxislabels = c("0", "0.25", "0.5", "0.75", "1"),  # smaller axis range
#           # seg = 4,
#           vlabels = paste0("  ", colnames(df)), # adds a small space to push labels outward
#           # paxislabels = rep("", length(colnames(df))),
#           # paxispos = rep(8.15, length(colnames(df))),
#           vlcex = 0.5,
#           title = title
#     )
#     par(xpd = TRUE)  # allow drawing outside plot area
#     legend(
#       x = "topright",               # or "bottomright", "topleft", etc.
#       # inset = c(-0.25, 0),          # move legend slightly outside plot area
#       legend = rownames(df),
#       bty = "n",
#       pch = 20,
#       col = colors_border[1:nrow(df)],
#       text.col = "black",
#       cex = 1.2,
#       pt.cex = 1.5
#     )
#     })
#   
#   # Return the plot as a grob object
#   return(p)
# }

# library(grid)
# library(gridExtra)
# library(png)
# 
# # helper to save a radar to a PNG and return a rasterGrob
# make_radar_png_grob <- function(df, title, filename = tempfile(fileext = ".png"),
#                                 width = 800, height = 800, res = 150) {
#   png(filename, width = width, height = height, res = res)
#   par(mar = c(1,1,3,1))
#   # draw the radarchart (use same code as before)
#   df_radar <- rbind(
#     max = rep(max(df, na.rm = TRUE) + 0.02, ncol(df)),
#     min = rep(0, ncol(df)),
#     df
#   )
#   df_radar_jitter <- df_radar
#   if(nrow(df_radar_jitter) > 2) {
#     df_radar_jitter[-c(1,2), ] <- df_radar_jitter[-c(1,2), ] +
#       matrix(runif((nrow(df_radar)-2) * ncol(df_radar), -0.02, 0.02),
#              ncol = ncol(df_radar))
#   }
#   
#   radarchart(
#     df_radar_jitter,
#     pcol = colors_border[1:nrow(df)],
#     plwd = 2,
#     plty = 1,
#     cglcol = "grey",
#     cglty = 1,
#     axislabcol = "grey",
#     vlabels = paste0("  ", colnames(df)),
#     vlcex = 0.5,
#     title = title
#   )
#   par(xpd = TRUE)
#   legend(
#     x = "topleft",
#     legend = rownames(df),
#     bty = "n",
#     pch = 20,
#     col = colors_border[1:nrow(df)],
#     text.col = "black",
#     cex = 0.6,
#     pt.cex = 0.8
#   )
#   dev.off()
#   
#   # read png back and convert to grob
#   img <- png::readPNG(filename)
#   rasterGrob(img, interpolate = TRUE)
# }

library(dplyr)
library(scales)
library(ggradar)
library(ggplot2)
library(tibble)

make_ggradar_plot <- function(df, 
                              title,
                              group_colours = colors_border,
                              axis_label_size = 2,
                              grid_label_size = 2.5,
                              legend_text_size = 12,       # smaller legend text
                              legend_key_spacing = 0.3) {
  
  # Ensure 'group' column exists
  if (!"group" %in% colnames(df)) {
    df <- df %>% tibble::rownames_to_column("group")
  }
  
  # Extract numeric columns and compute global min/max
  numeric_cols <- df %>% dplyr::select(where(is.numeric))
  global_min <- min(as.matrix(numeric_cols), na.rm = TRUE)
  global_max <- max(as.matrix(numeric_cols), na.rm = TRUE)
  
  # # Normalize all numeric columns using *global* min/max (not per column)
  df_norm <- df %>%
    dplyr::mutate(across(where(is.numeric),
                         ~ (. - global_min) / (global_max - global_min)))
  
  # Prepare grid labels in the original scale
  grid_breaks_original <- round(seq(global_min, global_max, length.out = 3), 2)
  names(grid_breaks_original) <- c("min", "mid", "max")
  
  # --- Create radar plot ---
  p <- ggradar::ggradar(
    df_norm,
    group.colours = group_colours,
    grid.min = 0,
    grid.mid = 0.5,
    grid.max = 1,
    values.radar = as.character(grid_breaks_original),
    background.circle.colour = "white",
    gridline.min.colour = "grey90",
    gridline.mid.colour = "grey80",
    gridline.max.colour = "grey70",
    axis.label.size = axis_label_size,
    grid.label.size = grid_label_size,
    group.point.size = 3,
    group.line.width = 1.2,
    label.gridline.min = TRUE,
    label.gridline.mid = TRUE,
    label.gridline.max = TRUE
  ) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),  # smaller legend text
      legend.key.size = unit(legend_key_spacing, "lines"),  # tighter spacing
      panel.background = element_blank()
    )
  
  return(p)
}

g1 <- make_ggradar_plot(mean_dcasted.list$AUC, title = "AUC")
g2 <- make_ggradar_plot(mean_dcasted.list$FeatNum, title = "Number of features")

###########

# heatmaps train/test features

###########
# load indics species result table for analysis_80_20
indics_species.df <- read.csv(file="~/Downloads/Mysubmissions/aime_07052025/data/projects/aime/analyses/bruvs/video_data_analysis/1.Article_figures/fig2/fig2-B-C-E-F/indics_species_analysis_results_80_20.csv")

# Helper function for generating Venn diagrams
generate_venn_plot <- function(venn_data, title_text, set_size = 3) {
  ggVennDiagram(venn_data, set_size = set_size, label_size = 5) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(color = "black", size = 9),
      legend.position = "bottom"
    ) + guides(fill = guide_legend(title = "IndicSpecies")) +
    scale_fill_distiller(palette = "RdBu", limits=c(0,50))
}

# filter indics species list by inshore_offshore at revalence 10%
indics_species_IN_OFF_prev_10.df <- indics_species.df[indics_species.df$prev_rate=="prev_10" & indics_species.df$comparison=="INSHORE_OFFSHORE",]

# rename bin and ter with predomics_ter and predomics_bin
indics_species_IN_OFF_prev_10.df[indics_species_IN_OFF_prev_10.df == 'terinter'] <- 'predomics_ter'
indics_species_IN_OFF_prev_10.df[indics_species_IN_OFF_prev_10.df == 'bininter'] <- 'predomics_bin'

library(ggVennDiagram)

# Precompute unique combinations
combos <- unique(indics_species_IN_OFF_prev_10.df[, c("seed", "prev_rate", "data", "source")])

# Initialize lists
venn_list_by_method <- list()
venn_list_by_source <- list()
venn_by_method <- list()
venn_by_source <- list()
input_data <- list()
common_indics_species <- list()

library(ggplot2)
# ---- Create venn_list_by_method and venn_list_by_source ----
for (r in seq_len(nrow(combos))) {
  i <- combos$seed[r]
  j <- combos$prev_rate[r]
  k <- combos$data[r]
  l <- combos$source[r]
  
  subset_df <- indics_species_IN_OFF_prev_10.df[
    indics_species_IN_OFF_prev_10.df$seed == i &
      indics_species_IN_OFF_prev_10.df$prev_rate == j &
      indics_species_IN_OFF_prev_10.df$data == k &
      indics_species_IN_OFF_prev_10.df$source == l,
  ]
  
  feats <- unique(subset_df$feature)
  
  venn_list_by_method[[i]][[j]][[k]][[l]] <- feats
  venn_list_by_source[[i]][[j]][[l]][[k]] <- feats
}

venn_by_method <- setNames(
  lapply(unique(combos$seed), function(i) {
    setNames(
      lapply(unique(combos$prev_rate), function(j) {
        setNames(
          lapply(unique(combos$data), function(k) {
            generate_venn_plot(
              venn_list_by_method[[i]][[j]][[k]],
              title_text = paste(i, j, k, sep = "/")
            )
          }),
          unique(combos$data) # names for k level
        )
      }),
      unique(combos$prev_rate) # names for j level
    )
  }),
  unique(combos$seed) # names for i level
)

# ---- Generate Venn plots by source ----
venn_by_source <- setNames(
  lapply(unique(combos$seed), function(i) {
    setNames(
      lapply(unique(combos$prev_rate), function(j) {
        setNames(
          lapply(unique(combos$source), function(k) {
            generate_venn_plot(
              venn_list_by_source[[i]][[j]][[k]],
              title_text = paste(i, j, k, sep = "/"),
              set_size = 2.5
            )
          }),
          unique(combos$source) # names for k level
        )
      }),
      unique(combos$prev_rate) # names for j level
    )
  }),
  unique(combos$seed) # names for i level
)

# ---- Input data extraction (simplified) ----
input_data <- setNames(
  lapply(unique(combos$seed), function(i) {
    setNames(
      lapply(unique(combos$prev_rate), function(j) {
        list(
          INSHORE_OFFSHORE = alldatalist.video_80_20[[i]][[j]][["terinter"]][["predout.maxn"]][["INSHORE_OFFSHORE"]][["Comp_data"]]
        )
      }),
      unique(combos$prev_rate)  # names for j level
    )
  }),
  unique(combos$seed)  # names for i level
)


# ---- Common indicator species ----
common_indics_species <- setNames(
  lapply(unique(combos$seed), function(i) {
    setNames(
      lapply(unique(combos$prev_rate), function(j) {
        setNames(
          lapply(unique(combos$data), function(k) {
            Reduce(intersect, venn_list_by_method[[i]][[j]][[k]])
          }),
          unique(combos$data)  # names for k level
        )
      }),
      unique(combos$prev_rate)  # names for j level
    )
  }),
  unique(combos$seed)  # names for i level
)


# Function to recursively flatten and bind nested lists
flatten_and_bind <- function(nested_list) {
  if (all(sapply(nested_list, is.data.frame))) {
    # If all elements are data frames, bind them
    do.call(rbind, nested_list)
  } else {
    # Otherwise, recursively process
    do.call(rbind, lapply(nested_list, flatten_and_bind))
  }
}
# heatmap for video data
load(file='~/Downloads/Mysubmissions/aime_07052025/data/projects/aime/analyses/bruvs/video_data_analysis/video_data_object.rda')

# video_data_abund <- sm$sample_info
# data_pres <- ifelse(sm$X>0, 1, 0)
# video_data_abund$richness <- colSums(sm$X>0)
sm$sample_info$Zone[sm$sample_info$hab %in% "BAR"] <- "Barrier"
sm$sample_info$Zone[sm$sample_info$hab %in% "BAY"] <- "Bay"
sm$sample_info$Zone[sm$sample_info$hab %in% "LAG"] <- "Lagoon"

#Get the raw data for the heatmap
dfaims <- sm$X
#Get the metadata for the heatmap
dfaims_meta <- sm$sample_info
#Get taxo info for species
dfaims_taxo <- sm$db_long
length(which(is.na(match(rownames(dfaims), dfaims_taxo$G_SP))))
dfaims_taxo <- dfaims_taxo[dfaims_taxo$G_SP %in% rownames(dfaims),]
dfaims_taxo <- unique(dfaims_taxo[,c("G_SP","Family","Genus","TrophicG")])
rownames(dfaims_taxo) <- dfaims_taxo$G_SP ; dfaims_taxo <- dfaims_taxo[,-1]
dfaims_sqrt <- apply(dfaims, 2, function(x){sqrt(x)})

### Make the heatmap on ggplot
dfaims_sqrt_melt <- as.data.frame(dfaims_sqrt)
dfaims_sqrt_melt$feature <- rownames(dfaims_sqrt_melt)
dfaims_sqrt_melt <- reshape2::melt(dfaims_sqrt_melt)
#Get the hclust ward.D on complete distance
dfaims_sqrt_clust.sp <- hclust(dist(dfaims_sqrt, method = "euclidean"), method = "ward.D") #cluster species with ward method from euclidean distances
dfaims_sqrt_clust.samples <- hclust(dist(t(dfaims_sqrt), method = "euclidean"), method = "ward.D") #cluster samples with ward method from euclidean distances
#Fix the samples and species order from the clustering results in the melted object
dfaims_sqrt_melt$feature <- factor(dfaims_sqrt_melt$feature, levels = dfaims_sqrt_clust.sp$labels[dfaims_sqrt_clust.sp$order])
dfaims_sqrt_melt$variable <- factor(dfaims_sqrt_melt$variable, levels = dfaims_sqrt_clust.samples$labels[dfaims_sqrt_clust.samples$order])
dfaims_sqrt_melt <- merge(dfaims_sqrt_melt, dfaims_meta[,c("hab","Habitat", "Zone")], by.x="variable", by.y=0, all.x=TRUE)
dfaims_taxo$feature <- rownames(dfaims_taxo)
dfaims_sqrt_melt <- merge(dfaims_sqrt_melt, dfaims_taxo[,c("Family","Genus", "TrophicG")], by.x="feature", by.y=0, all.x=TRUE)

pheat_map.plots <- list()
pheat_map_common.plots <- list()
common_dfaims_sqrt_melt.list <- list()
dfaims_sqrt_melt.list <- list()
for (i in unique(indics_species_IN_OFF_prev_10.df$data)) {
  for (j in unique(indics_species_IN_OFF_prev_10.df$prev_rate)) {
    for (k in unique(indics_species_IN_OFF_prev_10.df$seed)) {
      # Generate Venn diagram by source
      indics_species <- common_indics_species[[k]][[j]][[i]]
      training.X <- input_data[[k]][[j]][["INSHORE_OFFSHORE"]][["X_train"]]
      testing.X <- input_data[[k]][[j]][["INSHORE_OFFSHORE"]][["X_test"]]
      dfaims_sqrt_melt.list[[i]][[j]][[k]] <- dfaims_sqrt_melt
      dfaims_sqrt_melt.list[[i]][[j]][[k]]$prev_rate <- j
      dfaims_sqrt_melt.list[[i]][[j]][[k]]$data <- i
      dfaims_sqrt_melt.list[[i]][[j]][[k]]$seed <- k
      dfaims_sqrt_melt.list[[i]][[j]][[k]]$IsCommonIndSp <- ifelse(dfaims_sqrt_melt$feature %in% indics_species, 1, 0)
      dfaims_sqrt_melt.list[[i]][[j]][[k]]$sampling_set <- ifelse(dfaims_sqrt_melt$variable %in% rownames(training.X), "train",
                                                                  ifelse(dfaims_sqrt_melt$variable %in% rownames(testing.X), "test", NA))
      
      # get only common species
      common_dfaims_sqrt_melt.list[[i]][[j]][[k]] <- dfaims_sqrt_melt.list[[i]][[j]][[k]][dfaims_sqrt_melt.list[[i]][[j]][[k]]$feature %in% indics_species, ]
    }
    # bind all seed by rows
    #dfaims_sqrt_melt.list[[i]][[j]] <- rbind(dfaims_sqrt_melt.list[[i]][[j]][[k]])
  }
  dfaims_sqrt_melt.list[[i]] <- dfaims_sqrt_melt.df<- flatten_and_bind(dfaims_sqrt_melt.list[[i]])
  # reorder by seed
  dfaims_sqrt_melt.list[[i]] <- dfaims_sqrt_melt.list[[i]][order(dfaims_sqrt_melt.list[[i]]$seed, decreasing = FALSE), ]
  
  common_dfaims_sqrt_melt.list[[i]] <-  flatten_and_bind(common_dfaims_sqrt_melt.list[[i]])
  # reorder by seed
  common_dfaims_sqrt_melt.list[[i]] <- common_dfaims_sqrt_melt.list[[i]][order(common_dfaims_sqrt_melt.list[[i]]$seed, decreasing = FALSE), ]
  
  # plot pheatmap
  pheat_map.plots[[i]] <- ggplot(common_dfaims_sqrt_melt.list[[i]], aes(x=variable, y=feature, fill=value)) +
    geom_tile(colour=NA) +
    #scale_fill_viridis() +
    scale_fill_gradient(low = "gray", high = viridis::viridis(11)) +
    ylab("species") +
    xlab("samples") +
    labs(fill="Species Abundance") +
    ggtitle(paste("Inshore_Offshore", i, sep="/"))+
    # theme_minimal()+
    facet_grid(seed~Habitat+sampling_set, scales = "free_x", space = "free_x") +
    theme(strip.text = element_text(size = 10),
          strip.text.y = element_text(angle = 45, vjust = 0.5, hjust=0.7),
          axis.text.x = element_blank(),
          #axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  # plot pheatmap for common species identified by Indval and Predomics models
  pheat_map_common.plots[[i]] <- ggplot(common_dfaims_sqrt_melt.list[[i]], aes(x=variable, y=feature, fill=value)) +
    geom_tile(colour=NA) +
    #scale_fill_viridis() +
    scale_fill_gradient(low = "gray", high = viridis::viridis(11)) +
    ylab("species") +
    xlab("samples") +
    labs(fill="Species Abundance") +
    ggtitle(paste("Inshore_Offshore", i, sep="/"))+
    # theme_minimal()+
    facet_grid(seed~Habitat+sampling_set, scales = "free_x", space = "free_x") +
    theme(strip.text = element_text(size = 10),
          strip.text.y = element_text(angle = 45, vjust = 0.5, hjust=0.7),
          axis.text.x = element_blank(),
          #axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
}
library(ggh4x)
# pheatmap of indicator species by seed separated by train/test in Inshore-Offshore based on abundace
x <- pheat_map.plots[["maxN"]]$data
x.plots <- list()
for(i in unique(x$seed))
{
  print(i)
  x.plots[[i]] <- ggplot(x[x$seed %in% i,], aes(x=variable, y=feature, fill=value)) +
    geom_tile(colour=NA) +
    #scale_fill_viridis() +
    scale_fill_gradient(low = "white", high = viridis::viridis(11)) +
    ylab("species") +
    xlab("samples") +
    labs(fill="Species Abundance") +
    # ggtitle(paste("Inshore_Offshore", i, sep="/"))+
    # theme_minimal()+
    facet_nested(seed~Habitat+sampling_set, scales = "free_x", space = "free_x") +
    theme_bw() + 
    theme(strip.text = element_text(size = 10),
          strip.text.y = element_text(angle = 45),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

library(patchwork)
# create the heatmap combination
g3 <-wrap_plots(x.plots$seed_01 + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank()),
                x.plots$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                x.plots$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()), 
                x.plots$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()),
                ncol=1)
# save the updated version of Figure 5
pdf(file="/data/projects/aime/analyses/bruvs/Figure5_Update_EK/Figure5_AddRadar_plot.pdf", h=14, w=20)
hlay="
AACC
AACC
BBCC
BBCC
"
wrap_plots(wrap_elements(g1) + xlim(0, 100) + ylim(0, 50),
           wrap_elements(g2) + xlim(0, 100) + ylim(0, 50),
           wrap_elements(g3) + xlim(0, 100) + ylim(0, 50), 
           design = hlay) + plot_annotation(tag_levels = 'A')
dev.off()