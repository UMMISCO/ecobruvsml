library(predomics)
library(dplyr)
library(ggalt)
library(patchwork)
library(ggpubr)
library(ggVennDiagram)
library(reshape2)
library(data.table)
library(ggplot2)
library(ggradar)
library(fmsb)
library(gridExtra)

load(file = "analysis_scripts/all_analysis_results_80_20.Rda")

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
load(file = "analysis_scripts/rf_analysis_results_80_20.Rda")

# bind performances metrics for predomics best models and rf
metrics.list <- list(rf.metrics.df, pred_best_models.metrics.df)
metrics.df <- bind_rows(metrics.list)

# factorize methods
metrics.df$method <- factor(metrics.df$method, levels = c("BestModel_bininter", "BestModel_terinter", "RF_bininter","RF_Indval", "RF_terinter"))

# reorder seed from 01, 02, to 10
metrics.df= metrics.df[order(metrics.df$seed), ]

# remove rownames of metrics df
rownames(metrics.df) <- NULL

# #filter auc values from rf results
auc.df <- metrics.df[metrics.df$metrics=='auc', ]

#filter acc values from rf results
acc.df <- metrics.df[metrics.df$metrics=='accuracy', ]

#filter f1 values from rf results
f1.df <- metrics.df[metrics.df$metrics=='f1', ]

# select only perf for prev_10
auc_prev_10.df <- auc.df[auc.df$prev_rate=='prev_10', ]

# plot auc perf for prev_10
ggplot(data=auc_prev_10.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  geom_point() + 
  # geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "AUC") +
  facet_grid(.~source) +
  # ylim(c(0.6,1.05))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = after_stat(p.signif)), method="kruskal", label.y = c(1.01, 1.02, 1.03,1.04)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12), strip.text.x = element_text(size=12))+
  # scale_size(range = c(0.3,1)) +
  scale_y_continuous(breaks = seq(0.3, 1, 0.1))


auc_prev_10.dly.df <- plyr::ddply(auc_prev_10.df,
                            c("metrics","comparison","source","prev_rate","method"), 
                            summarise,
                            mean.value=mean(value), 
                            # val.na=length(which(is.na(value))),
                            sd.value=sd(value), 
                            se.value=sd(value)/sqrt(length(value)),
                            mean.feats=mean(featNum),
                            sd.feats=sd(featNum),
                            se.feats=sd(featNum)/sqrt(length(featNum)))

library(ggh4x)
p1 <- ggplot(data=auc_prev_10.dly.df, aes(x= method, y= mean.value))+
  # geom_boxplot(alpha = 0.6)+
  geom_point() + 
  ylab("AUC\n(mean \u00B1 std.err)") + 
  geom_errorbar(aes(ymin=mean.value-se.value, ymax=mean.value+se.value), width = .2) + 
  facet_nested(.~source+comparison) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=6.5))


p2 <- ggplot(data=auc_prev_10.dly.df, aes(x= method, y= mean.feats))+
  # geom_boxplot(alpha = 0.6)+
  geom_point() + 
  geom_errorbar(aes(ymin=mean.feats-se.feats, ymax=mean.feats+se.feats), width = .2) + 
  facet_nested(.~source+comparison) + 
  ylab("Features\n(mean \u00B1 std.err)") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=6.5))

##pairwise comparisons of AUC values btw methods
dunnlist <- list()
library(dunn.test)
for(i in unique(auc_prev_10.df$comparison))
{
  for(j in unique(auc_prev_10.df$source))
  {
    ijdf <- auc_prev_10.df[auc_prev_10.df$comparison %in% i & auc_prev_10.df$source %in% j,]
    ijdf.dunn <- as.data.frame(dunn.test(x = ijdf$value, g = ijdf$method, kw = TRUE, method = "none", table = TRUE))
    ijdf.dunn$comparison <- i
    ijdf.dunn$source <- j
    dunnlist[[paste(i,j)]] <- ijdf.dunn
  }
}
dunnlist.df <- do.call("rbind", dunnlist)
##Plot dumbell from pairs contrast
dunnlist.df$comparisons_A <- gsub(" .*","", dunnlist.df$comparisons)
dunnlist.df$comparisons_B <- gsub("^.* ","", dunnlist.df$comparisons)
dunnlist.df$comparisons_A <- as.numeric(as.character(factor(dunnlist.df$comparisons_A, levels = unique(auc_prev_10.df$method)[order(unique(auc_prev_10.df$method))], labels = seq(1:length(unique(auc_prev_10.df$method))))))
dunnlist.df$comparisons_B <- as.numeric(as.character(factor(dunnlist.df$comparisons_B, levels = unique(auc_prev_10.df$method)[order(unique(auc_prev_10.df$method))], labels = seq(1:length(unique(auc_prev_10.df$method))))))
dunnlist.df$P.adjusted.cat <- ifelse(dunnlist.df$P.adjusted<0.05,"FDR<0.05","FDR>0.05")

p1_dumbell <- ggplot(dunnlist.df[dunnlist.df$P.adjusted<0.05,], aes(x = comparisons_A, xend = comparisons_B, y = comparisons,group = comparisons)) +
  geom_dumbbell(aes(colour=P.adjusted.cat)) + 
  scale_x_continuous(breaks = seq(1,length(unique(auc_prev_10.df$method))), labels = unique(auc_prev_10.df$method)[order(unique(auc_prev_10.df$method))]) + 
  facet_nested(.~source+comparison) + 
  scale_colour_manual(values = c("black","grey")) + 
  xlab("Models") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(fill = NA, linetype = 0),
        legend.position = "none")

hlay <- "
A
A
B
"
p1 + p1_dumbell + plot_layout(design = hlay)

##pairwise comparison model sizes btw models 
dunnlist2 <- list()
library(dunn.test)
for(i in unique(auc_prev_10.df$comparison))
{
  for(j in unique(auc_prev_10.df$source))
  {
    ijdf <- auc_prev_10.df[auc_prev_10.df$comparison %in% i & auc_prev_10.df$source %in% j,]
    ijdf.dunn <- as.data.frame(dunn.test(x = ijdf$featNum, g = ijdf$method, kw = TRUE, method = "none", table = TRUE))
    ijdf.dunn$comparison <- i
    ijdf.dunn$source <- j
    dunnlist2[[paste(i,j)]] <- ijdf.dunn
  }
}
dunnlist2.df <- do.call("rbind", dunnlist2)
##Plot dumbell from pairs contrast
dunnlist2.df$comparisons_A <- gsub(" .*","", dunnlist2.df$comparisons)
dunnlist2.df$comparisons_B <- gsub("^.* ","", dunnlist2.df$comparisons)
dunnlist2.df$comparisons_A <- as.numeric(as.character(factor(dunnlist2.df$comparisons_A, levels = unique(auc_prev_10.df$method)[order(unique(auc_prev_10.df$method))], labels = seq(1:length(unique(auc_prev_10.df$method))))))
dunnlist2.df$comparisons_B <- as.numeric(as.character(factor(dunnlist2.df$comparisons_B, levels = unique(auc_prev_10.df$method)[order(unique(auc_prev_10.df$method))], labels = seq(1:length(unique(auc_prev_10.df$method))))))
dunnlist2.df$P.adjusted.cat <- ifelse(dunnlist2.df$P.adjusted<0.05,"FDR<0.05","FDR>0.05")

p2_dumbell <- ggplot(dunnlist2.df[dunnlist2.df$P.adjusted<0.05,], aes(x = comparisons_A, xend = comparisons_B, y = comparisons,group = comparisons)) +
  geom_dumbbell(aes(colour=P.adjusted.cat)) + 
  scale_x_continuous(breaks = seq(1,length(unique(auc_prev_10.df$method))), labels = unique(auc_prev_10.df$method)[order(unique(auc_prev_10.df$method))]) + 
  facet_nested(.~source+comparison) + 
  xlab("Models") + 
  scale_colour_manual(values = c("black","grey")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(fill = NA, linetype = 0),
        legend.position = "none")

###########

# heatmaps train/test features

###########
# load result list of 80_20 analysis (Indval/terinter/bininter)
load("analysis_scripts/all_analysis_results_80_20.Rda")

names(alldatalist.video_80_20) <- ifelse(names(alldatalist.video_80_20)=="seed_10", "seed_10", gsub("seed_", "seed_0", names(alldatalist.video_80_20)))

# load indics species result table for analysis_80_20
indics_species.df <- read.csv(file="analysis_scripts/indics_species_analysis_results_80_20.csv")

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

# # filter indics species list by inshore_offshore
# indics_species.df <- indics_species.df[indics_species.df$comparison=="INSHORE_OFFSHORE",]

# rename bin and ter with predomics_ter and predomics_bin
indics_species.df[indics_species.df == 'terinter'] <- 'predomics_ter'
indics_species.df[indics_species.df == 'bininter'] <- 'predomics_bin'

# get only species at revalence 10%

indics_species_prev_10.df <- indics_species.df[indics_species.df$prev_rate=="prev_10",]

# Get unique values for each dimension
seeds <- unique(indics_species_prev_10.df$seed)
prev_rates <- unique(indics_species_prev_10.df$prev_rate)
data_types <- unique(indics_species_prev_10.df$data)
sources <- unique(indics_species_prev_10.df$source)
comparisons <- unique(indics_species_prev_10.df$comparison)

# Function to get features for specific combination
get_features <- function(seed_val, prev_val, data_val, source_val, comparison_val) {
  unique(indics_species_prev_10.df$feature[
    (indics_species_prev_10.df$seed == seed_val) &
      (indics_species_prev_10.df$prev_rate == prev_val) &
      (indics_species_prev_10.df$data == data_val) &
      (indics_species_prev_10.df$source == source_val) &
      (indics_species_prev_10.df$comparison == comparison_val)
  ])
}

# Venn lists by method

venn_list_by_method <- lapply(seeds, function(i) {
  seed_list <- lapply(prev_rates, function(j) {
    prev_list <- lapply(data_types, function(k) {
      data_list <- lapply(sources, function(l) {
        comp_list <- lapply(comparisons, function(m) {
          get_features(i, j, k, l, m)
        })
        names(comp_list) <- comparisons
        comp_list
      })
      names(data_list) <- sources
      data_list
    })
    names(prev_list) <- data_types
    prev_list
  })
  names(seed_list) <- prev_rates
  seed_list
})
names(venn_list_by_method) <- seeds

# Venn lists by source

venn_list_by_source <- lapply(seeds, function(i) {
  seed_list <- lapply(prev_rates, function(j) {
    source_list <- lapply(sources, function(l) {
      data_list <- lapply(data_types, function(k) {
        venn_list_by_method[[i]][[j]][[k]][[l]]
      })
      names(data_list) <- data_types
      data_list
    })
    names(source_list) <- sources
    source_list
  })
  names(seed_list) <- prev_rates
  seed_list
})
names(venn_list_by_source) <- seeds

# Generate Venn diagrams by method

venn_by_method <- lapply(seeds, function(i) {
  seed_list <- lapply(prev_rates, function(j) {
    data_list <- lapply(data_types, function(k) {
      comp_list <- lapply(comparisons, function(m) {
        
        # FIXED: must pass a list of sets to ggVennDiagram
        sets_list <- lapply(sources, function(l) {
          venn_list_by_method[[i]][[j]][[k]][[l]][[m]]
        })
        names(sets_list) <- sources
        
        generate_venn_plot(
          sets_list,
          title_text = paste(i, j, k, m, sep = "/")
        )
      })
      names(comp_list) <- comparisons
      comp_list
    })
    names(data_list) <- data_types
    data_list
  })
  names(seed_list) <- prev_rates
  seed_list
})
names(venn_by_method) <- seeds

# Generate Venn diagrams by source

venn_by_source <- lapply(seeds, function(i) {
  seed_list <- lapply(prev_rates, function(j) {
    source_list <- lapply(sources, function(l) {
      comp_list <- lapply(comparisons, function(m) {
        
        sets_list <- lapply(data_types, function(k) {
          venn_list_by_source[[i]][[j]][[l]][[k]][[m]]
        })
        names(sets_list) <- data_types
        
        generate_venn_plot(
          sets_list,
          title_text = paste(i, j, l, m, sep = "/"),
          set_size = 2.5
        )
      })
      names(comp_list) <- comparisons
      comp_list
    })
    names(source_list) <- sources
    source_list
  })
  names(seed_list) <- prev_rates
  seed_list
})
names(venn_by_source) <- seeds

# Input data by comparison

input_data <- lapply(seeds, function(i) {
  seed_list <- lapply(prev_rates, function(j) {
    comp_list <- lapply(comparisons, function(m) {
      alldatalist.video_80_20[[i]][[j]][["terinter"]][["predout.maxn"]][[m]][["Comp_data"]]
    })
    names(comp_list) <- comparisons
    comp_list
  })
  names(seed_list) <- prev_rates
  seed_list
})
names(input_data) <- seeds

# Common indicator species

common_indics_species <- lapply(prev_rates, function(j) {
  prev_list <- lapply(data_types, function(k) {
    data_list <- lapply(comparisons, function(l) {
      
      # For each source, first intersect across all seeds
      source_intersections <- lapply(sources, function(m) {
        
        # Species list for each seed
        species_by_seed <- lapply(seeds, function(i) {
          unlist(venn_list_by_method[[i]][[j]][[k]][[m]][[l]])
        })
        
        # Union across seeds for ONE source
        Reduce(union, species_by_seed)
      })
      names(source_intersections) <- sources
      
      # Now intersect across sources
      all_sources_intersection <- Reduce(intersect, source_intersections)
      
      all_sources_intersection
    })
    names(data_list) <- comparisons
    data_list
  })
  names(prev_list) <- data_types
  prev_list
})
names(common_indics_species) <- prev_rates

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
load(file='data/video_data_object.rda')

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
dfaims_sqrt_melt <- melt(dfaims_sqrt_melt)
#Get the hclust ward.D on complete distance
dfaims_sqrt_clust.sp <- hclust(dist(dfaims_sqrt, method = "euclidean"), method = "ward.D") #cluster species with ward method from euclidean distances
dfaims_sqrt_clust.samples <- hclust(dist(t(dfaims_sqrt), method = "euclidean"), method = "ward.D") #cluster samples with ward method from euclidean distances
#Fix the samples and species order from the clustering results in the melted object
dfaims_sqrt_melt$feature <- factor(dfaims_sqrt_melt$feature, levels = dfaims_sqrt_clust.sp$labels[dfaims_sqrt_clust.sp$order])
dfaims_sqrt_melt$variable <- factor(dfaims_sqrt_melt$variable, levels = dfaims_sqrt_clust.samples$labels[dfaims_sqrt_clust.samples$order])
dfaims_sqrt_melt <- merge(dfaims_sqrt_melt, dfaims_meta[,c("hab","Habitat", "Zone")], by.x="variable", by.y=0, all.x=TRUE)
dfaims_taxo$feature <- rownames(dfaims_taxo)
dfaims_sqrt_melt <- merge(dfaims_sqrt_melt, dfaims_taxo[,c("Family","Genus", "TrophicG")], by.x="feature", by.y=0, all.x=TRUE)

# Initialize top-level lists
dfaims_sqrt_melt.list <- list()
common_dfaims_sqrt_melt.list <- list()

# Loop over data types
dfaims_sqrt_melt.list <- lapply(data_types, function(i) {
  
  prev_list <- lapply(prev_rates, function(j) {
    
    seed_list <- lapply(seeds, function(k) {
      
      # Loop over comparisons
      comp_list <- lapply(comparisons, function(m) {
        
        # common indicator species
        indics_species <- common_indics_species[[j]][[i]][[m]]
        
        # train/test sets for this comparison
        training.X <- input_data[[k]][[j]][[m]][["X_train"]]
        testing.X <- input_data[[k]][[j]][[m]][["X_test"]]
        
        # copy dfaims_sqrt_melt and add metadata
        df_copy <- dfaims_sqrt_melt
        df_copy$prev_rate <- j
        df_copy$data <- i
        df_copy$seed <- k
        df_copy$comparison <- m
        df_copy$IsCommonIndSp <- ifelse(df_copy$feature %in% indics_species, 1, 0)
        df_copy$sampling_set <- ifelse(df_copy$variable %in% rownames(training.X), "train",
                                       ifelse(df_copy$variable %in% rownames(testing.X), "test", NA))
        
        # Remove rows that do not belong to train or test
        df_copy <- df_copy[!is.na(df_copy$sampling_set), ]
        
        df_copy
      })
      
      names(comp_list) <- comparisons
      comp_list
    })
    
    names(seed_list) <- seeds
    seed_list
  })
  
  names(prev_list) <- prev_rates
  prev_list
})

names(dfaims_sqrt_melt.list) <- data_types

# -----------------------------
# Build common species list
# -----------------------------
common_dfaims_sqrt_melt.list <- lapply(dfaims_sqrt_melt.list, function(prev_list) {
  lapply(prev_list, function(seed_list) {
    lapply(seed_list, function(comp_list) {
      lapply(comp_list, function(df) df[df$IsCommonIndSp == 1, ])
    })
  })
})

# -----------------------------
# Flatten by data type and order by seed
# -----------------------------
dfaims_sqrt_melt.list <- lapply(dfaims_sqrt_melt.list, function(prev_list) {
  flat <- flatten_and_bind(prev_list)
  flat[order(flat$seed), ]
})

common_dfaims_sqrt_melt.list <- lapply(common_dfaims_sqrt_melt.list, function(prev_list) {
  flat <- flatten_and_bind(prev_list)
  flat[order(flat$seed), ]
})


# Extract the data for maxN
x <- common_dfaims_sqrt_melt.list[["maxN"]]

# Generate heatmaps with lapply (no loops)
# Unique values
all_comparisons <- unique(x$comparison)
all_seeds       <- unique(x$seed)

x.plots <- setNames(
  lapply(all_comparisons, function(comp) {
    
    comp_subset <- x[x$comparison == comp, ]
    
    # Choose Habitat column depending on comparison
    facet_var <- if (comp == "INSHORE_OFFSHORE") "Habitat" else "hab"
    
    setNames(
      lapply(all_seeds, function(seed_val) {
        
        df_subset <- comp_subset[comp_subset$seed == seed_val, ]
        
        ggplot(df_subset, aes(x = variable, y = feature, fill = value)) +
          geom_tile(colour = NA) +
          scale_fill_gradient(low = "white", high = viridis::viridis(11)) +
          ylab("species") +
          xlab("samples") +
          labs(fill = "Species Abundance") +
          facet_nested(
            seed ~ .data[[facet_var]] + sampling_set,
            scales = "free_x",
            space = "free_x"
          ) +
          theme_bw() +
          theme(
            strip.text = element_text(size = 10),
            strip.text.y = element_text(angle = 360),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
      }),
      all_seeds
    )
  }),
  all_comparisons
)

hlay1="
A
A
B
"
g1 <- wrap_plots(p1,p1_dumbell+ theme(axis.text.y = element_blank()), design = hlay1)
g2 <- wrap_plots(p2,p2_dumbell+ theme(axis.text.y = element_blank()), design = hlay1)
g3 <-wrap_plots(x.plots$INSHORE_OFFSHORE$seed_01 + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                x.plots$INSHORE_OFFSHORE$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                x.plots$INSHORE_OFFSHORE$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$INSHORE_OFFSHORE$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                ncol=1)
# save the Figure 5
pdf(file="figures/Figure5/Figure5_old.pdf", h=8, w=16)
hlay="
AAAACC
AAAACC
BBBBCC
BBBBCC
"
wrap_plots(wrap_elements(g1) + xlim(0, 100) + ylim(0, 50),
           wrap_elements(g2) + xlim(0, 100) + ylim(0, 50),
           wrap_elements(g3) + xlim(0, 100) + ylim(0, 50), 
           design = hlay) + plot_annotation(tag_levels = 'A')
dev.off()

# save the Figure S6
# pdf(file="figures/Figure5/FigureS6.pdf", h=12, w=12)

pdf(file="figures/Figure5/FigureS6.pdf", h=12, w=12)
hlay2="
A
A
B
C
C
D
"
p1 + ggtitle("A") + 
  p1_dumbell + theme(axis.text.y = element_blank()) +
  p2 + ggtitle("B") + 
  p2_dumbell + theme(axis.text.y = element_blank()) + 
  plot_layout(design = hlay2)

dev.off()

comp1 <- wrap_plots(x.plots$BAR_BAY$seed_01 + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    x.plots$BAR_BAY$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    x.plots$BAR_BAY$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_BAY$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    ncol=1)

comp2 <- wrap_plots(x.plots$BAR_LAG$seed_01 + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    x.plots$BAR_LAG$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    x.plots$BAR_LAG$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAR_LAG$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    ncol=1)

comp3 <- wrap_plots(x.plots$BAY_LAG$seed_01 + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    x.plots$BAY_LAG$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    x.plots$BAY_LAG$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                    x.plots$BAY_LAG$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                    ncol=1)

pdf(file="figures/Figure5/FigureS4.pdf", h=8, w=16)
hlay3=
  "
AABBCC
AABBCC
AABBCC
"
wrap_plots(wrap_elements(comp1),
           wrap_elements(comp2),
           wrap_elements(comp3), 
           design = hlay3) + plot_annotation(tag_levels = 'A')
dev.off()


###########

# Radar plots of AUC performances of models

###########

# filtering AUC metrics 
metrics_filtered.df <- metrics.df[metrics.df$metrics=="auc" & metrics.df$prev_rate == "prev_10",]

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

# Define colors
colors_border <- c("#9467bd", "#ff7f0e", "#648C67", "#8195B2", "#CC3837")
colors_in <- scales::alpha(colors_border, 0.15)

make_ggradar_plot <- function(df, 
                              title,
                              group_colours = colors_border,
                              axis_label_size = 2,
                              grid_label_size = 2.5,
                              legend_text_size = 12,       
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

p1 <- make_ggradar_plot(mean_dcasted.list$AUC, title = "AUC")
p2 <- make_ggradar_plot(mean_dcasted.list$FeatNum, title = "Number of features")

# save the updated version of Figure 5 with radar plot
pdf(file="figures/Figure5/Figure5_AddRadar_plot.pdf", h=14, w=20)
hlay4="
AACC
AACC
BBCC
BBCC
"
wrap_plots(wrap_elements(p1) + xlim(0, 100) + ylim(0, 50),
           wrap_elements(p2) + xlim(0, 100) + ylim(0, 50),
           wrap_elements(g3) + xlim(0, 100) + ylim(0, 50), 
           design = hlay4) + plot_annotation(tag_levels = 'A')
dev.off()