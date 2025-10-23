library(predomics)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggh4x)
library(ggVennDiagram)
library(reshape2)
library(patchwork)
library(tidyr)


load(file = "/data/projects/aime/analyses/bruvs/FigureS5/results_replicas_analysis/all_analysis_results_split_replicas.Rda")

# Define parameters
params <- expand.grid(
  i = names(alldatalist.video_split_replicas),
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
  resultObj <- alldatalist.video_split_replicas[[i]][[j]][[k]][[l]][[m]]
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
  
  # # Safely round metrics
  # safe_round <- function(x, digits = 3) {
  #   if (is.null(x) || length(x) == 0 || !is.numeric(x)) return(NA_real_)
  #   round(as.numeric(x), digits)
  # }
  # 
  # # Get metrics (AUC, ACC, f1)
  # best.model.auc <- safe_round(best.model.test$auc_)
  # best.model.acc <- safe_round(best.model.test$accuracy_)
  # best.model.f1 <- safe_round(best.model.test$f1_)
  
  
  # Store metrics
  data.frame(
    metrics = c("auc", "accuracy", "f1"),
    value = c(best.model.auc, best.model.acc, best.model.f1),
    comparison = m,
    source = ifelse(l == "predout.bin", "presAbs", "maxN"),
    featNum = length(best.model$names_),
    prev_rate = j,
    method = paste("BestModel", k, sep = "_"),
    replica_comp= i
  )
}

# Apply the function to each row of parameter combinations
results <- lapply(seq_len(nrow(params)), function(row_idx) {
  row <- params[row_idx, ]
  best_model_avaluation(row)
})

# bind list into a dafarame
pred_best_models.metrics.df <- bind_rows(results)

# load rf analysis results on split samples by replicas for Indval and Predomics Models
load(file = "/data/projects/aime/analyses/bruvs/FigureS5/rf_analysis_by_replicas/rf_analysis_results_split_by_replicas.Rda")

# load rf analysis results on split samples by replicas for Predomics Best Models
load(file = "/data/projects/aime/analyses/bruvs/FigureS5/rf_analysis_by_replicas/rf_analysis_results_split_by_replicas_BestModel_Predomics.Rda")

# bind performances metrics for predomics best models and rf
#metrics.list <- list(rf.metrics.df, pred_best_models.metrics.df)
metrics.df <- bind_rows(rf.metrics.df, pred_best_models.metrics.df, rf.metrics.bestModel.df)

# factorize methods
metrics.df$method <- factor(metrics.df$method, levels = c("BestModel_bininter", "BestModel_terinter", "RF_bestModel_bininter", "RF_bestModel_terinter", "RF_bininter","RF_Indval", "RF_terinter"))

# factorise replica_comp
metrics.df$replica_comp <- ifelse(metrics.df$replica_comp=="ABORE_vs_MBERE", "train on ABORE and\ntest on MBERE", "train on MBERE and\ntest on ABORE")

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

p1 <- ggplot(data=auc_prev_10.df, aes(x=method, y=value, fill = featNum)) +
  geom_bar(stat="identity", position = "dodge")+
  # ylab("AUC") +
  # xlab("Models") +
  labs(x="Models", y= "AUC", fill= "Number of species") +
  scale_fill_viridis(option = "D") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_text(angle= 90, margin = ggplot2::margin(r = 15, b=10)),
        strip.text.y = element_text(angle = 90),
        # strip.text.x = element_text(size=8)
  )+
  #scale_fill_gradient(low = "black", high = "grey") +
  # scale_x_discrete(limits = c("RF_Indval", "RF_terinter", "RF_bininter", "RF_bestModel_terinter", "RF_bestModel_bininter", "BestModel_terinter", "BestModel_bininter")) +
  facet_nested(replica_comp~source+comparison)

###########

# heatmaps train/test features

###########

# load result list of replicas analysis for Best Models (Indval/terinter/bininter)
load(file = "/data/projects/aime/analyses/bruvs/FigureS5/results_replicas_analysis/all_analysis_results_split_replicas_BestModels_Predomics.Rda")

# load indics species result table for replicas analysis with bininter/indval/terinter
indics_species.df <- read.csv(file="/data/projects/aime/analyses/bruvs/FigureS5/results_replicas_analysis/indics_species_analysis_split_by_replicas.csv")

# load indics species result table for replicas analysis with Best models of bininter/terinter
indics_species_BM.df <- read.csv(file="/data/projects/aime/analyses/bruvs/FigureS5/results_replicas_analysis/indics_species_analysis_split_by_replicas_BestModels_Predomics.csv")

# get indics species for all models
all_indics_species.df <- rbind(indics_species.df, indics_species_BM.df)
# all_indics_species_prev_10.df <- all_indics_species.df[all_indics_species.df$prev_rate %in% "prev_10",]

# Initialize lists
venn_list_by_method <- list()
venn_by_method <- list()
venn_list_by_source <- list()
venn_by_source <- list()

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

# filter indics species list by inshore_offshore
# indics_species_IN_OFF.df <- all_indics_species.df[all_indics_species.df$comparison=="INSHORE_OFFSHORE",]

# rename bin and ter with predomics_ter and predomics_bin
all_indics_species.df[all_indics_species.df == 'terinter'] <- 'Predomics_ter'
all_indics_species.df[all_indics_species.df == 'bininter'] <- 'Predomics_bin'
all_indics_species.df[all_indics_species.df == 'BestModel_terinter'] <- 'BestModel_Predomics_ter'
all_indics_species.df[all_indics_species.df == 'BestModel_bininter'] <- 'BestModel_Predomics_bin'

# get only species at Prevalence 10%

# indics_species_IN_OFF_prev_10.df <- all_indics_species.df[all_indics_species.df$prev_rate=="prev_10",]

#Main loops to generate the Venn diagrams
for (i in unique(all_indics_species.df$replica_comp)) {
  for (j in unique(all_indics_species.df$comparison)) {
    for (k in unique(all_indics_species.df$prev_rate)) {
      for (l in unique(all_indics_species.df$data)) {
        for (m in unique(all_indics_species.df$source)) {
          # Create Venn list for method
          venn_list_by_method[[i]][[j]][[k]][[l]][[m]] <- unique(all_indics_species.df$feature[
            (all_indics_species.df$replica_comp == i) &
              (all_indics_species.df$comparison == j) &
              (all_indics_species.df$prev_rate == k) &
              (all_indics_species.df$data == l) &
              (all_indics_species.df$source == m)
          ])
          # Create Venn list for method
          venn_list_by_source[[i]][[j]][[k]][[m]][[l]] <- venn_list_by_method[[i]][[j]][[k]][[l]][[m]]
        }
        # Generate Venn diagram by method
        venn_by_method[[i]][[j]][[k]][[l]] <- generate_venn_plot(venn_list_by_method[[i]][[j]][[k]][[l]], title_text = paste(i,j,k, sep="/"))
      }
    }
  }
  
}

for (i in unique(all_indics_species.df$replica_comp)) {
  for (j in unique(all_indics_species.df$comparison)) {
    for (k in unique(all_indics_species.df$prev_rate)) {
      for (l in unique(all_indics_species.df$source)) {
        # Generate Venn diagram by source
        venn_by_source[[i]][[j]][[k]][[l]] <- generate_venn_plot(venn_list_by_source[[i]][[j]][[k]][[l]], title_text = paste(i,j,k,l, sep="/"), set_size = 2.5)
      }
    }
  }
}

input_data= list()

# get input data for each replica_comp for prev_0 and prev_10
for (i in unique(all_indics_species.df$replica_comp)) {
  for (j in unique(all_indics_species.df$prev_rate)) {
    for (k in unique(all_indics_species.df$comparison)) {
      for (l in unique(all_indics_species.df$data)) {
        if(l=="maxN"){
          input_data[[i]][[j]][[k]][[l]]= alldatalist.video_split_replicas[[i]][[j]][["terinter"]][["predout.maxn"]][[k]][["Comp_data"]]
        }else{
          input_data[[i]][[j]][[k]][[l]]= alldatalist.video_split_replicas[[i]][[j]][["terinter"]][["predout.bin"]][[k]][["Comp_data"]]
        }
        
      }
      
    }
  }
}

all_methods_indics_species <- list()
for (i in unique(all_indics_species.df$replica_comp)) {
  for (j in unique(all_indics_species.df$comparison)) {
    for (k in unique(all_indics_species.df$prev_rate)) {
      for (l in unique(all_indics_species.df$data)) {
        # Generate Venn diagram by source
        all_methods_indics_species[[l]][[k]][[j]][[i]] <- Reduce(union, venn_list_by_method[[i]][[j]][[k]][[l]])
      }
    }
  }
}


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
load(file='/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/video_data_object.rda')

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

flatten_recursive <- function(x) {
  if (is.list(x)) {
    return(unlist(lapply(x, flatten_recursive), recursive = FALSE))
  } else {
    return(list(x))
  }
}

pheat_map.plots <- list()
pheat_map_all.plots <- list()
all_dfaims_sqrt_melt.list <- list()
dfaims_sqrt_melt.list <- list()

for (i in unique(all_indics_species.df$data)) {
  for (j in unique(all_indics_species.df$prev_rate)) {
    for (k in unique(all_indics_species.df$comparison)) {
      
      # get the list of all indics species by data, prev_rate, comparison
      indics_species <- unique(unlist(flatten_recursive(all_methods_indics_species[[i]][[j]][[k]])))
      
      for (l in unique(all_indics_species.df$replica_comp)) {
        # indics_species <- unique(all_methods_indics_species[[i]][[j]][[l]][[k]])
        training.X <- input_data[[l]][[j]][[k]][[i]][["X_train"]]
        testing.X <- input_data[[l]][[j]][[k]][[i]][["X_test"]]
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]] <- dfaims_sqrt_melt
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$data <- i
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$prev_rate <- j
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$comparison <- k
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$replica_comp <- l
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$sampling_set <- ifelse(dfaims_sqrt_melt$variable %in% rownames(training.X), gsub("_.*", "", l),
                                                                         ifelse(dfaims_sqrt_melt$variable %in% rownames(testing.X), gsub(".*_", "", l), NA))
        # Remove NA if they exist from sampling set
        dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]] <- dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]][!is.na(dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$sampling_set), ]
        
        # get only common species
        all_dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]] <- dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]][dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]]$feature %in% indics_species, ]
      }
      
      # bind all seed by rows
      
      dfaims_sqrt_melt.list[[i]][[j]][[k]] <- rbind(dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]])
      
      all_dfaims_sqrt_melt.list[[i]][[j]][[k]] <-  rbind(all_dfaims_sqrt_melt.list[[i]][[j]][[k]][[l]])
      
      # plot pheatmap
      pheat_map.plots[[i]][[j]][[k]] <- ggplot(dfaims_sqrt_melt.list[[i]][[j]][[k]], aes(x=variable, y=feature, fill=value)) +
        geom_tile(colour=NA) +
        # scale_fill_viridis() +
        scale_fill_gradient(low = "white", high = viridis::viridis(11)) +
        ylab("species") +
        xlab("samples") +
        labs(fill="Species Abundance") +
        # ggtitle(paste("Inshore_Offshore", i, sep="/"))+
        # theme_minimal()+
        facet_nested(.~Zone+sampling_set, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme( #strip.text = element_text(size = 10),
          strip.text.y = element_text(angle = 45, vjust = 0.5, hjust=0.7),
          axis.text.x = element_blank(),
          # axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.text.y = element_blank(),
          # axis.ticks.y = element_blank(),
          legend.title = element_text(angle= 90, margin = ggplot2::margin(r = 25, b=10)),
          # legend.position = "none"
        )
      
      # plot pheatmap for common species identified by Indval and Predomics models
      pheat_map_all.plots[[i]][[j]][[k]] <- ggplot(all_dfaims_sqrt_melt.list[[i]][[j]][[k]], aes(x=variable, y=feature, fill=value)) +
        geom_tile(colour=NA) +
        #scale_fill_viridis() +
        scale_fill_gradient(low = "white", high = viridis::viridis(11)) +
        ylab("species") +
        xlab("samples") +
        labs(fill="Species Abundance") +
        # ggtitle(paste("Inshore_Offshore", i, sep="/"))+
        # theme_minimal()+
        facet_nested(.~Zone+sampling_set, scales = "free_x", space = "free_x") +
        theme_bw() +
        theme(# strip.text = element_text(size = 10),
          strip.text.y = element_text(angle = 45, vjust = 0.5, hjust=0.7),
          axis.text.x = element_blank(),
          # axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.text.y = element_blank(),
          # axis.ticks.y = element_blank(),
          legend.title = element_text(angle= 90, margin = ggplot2::margin(r = 25, b=10)),
          # legend.position = "none"
        )
    }
  }
}

pdf(file="/data/projects/aime/analyses/bruvs/FigureS5/FigureS5_Best.pdf", h=13, w=11)
hlay="
A
A
B
B
B
"
p1 + ggtitle("A") + 
  pheat_map_all.plots$maxN$prev_10$INSHORE_OFFSHORE + ggtitle("B") + 
  plot_layout(design = hlay)
dev.off()

pdf(file="/data/projects/aime/analyses/bruvs/FigureS5/FigureS5_Best_update.pdf", h=22, w=15)
hlay="
AAAA
BBCC
BBCC
DDEE
DDEE
"
p1 + ggtitle("A") + 
  pheat_map_all.plots$maxN$prev_10$INSHORE_OFFSHORE + ggtitle("B") + 
  pheat_map_all.plots$maxN$prev_10$BAR_BAY + ggtitle("C") + 
  pheat_map_all.plots$maxN$prev_10$BAR_LAG + ggtitle("D") + 
  pheat_map_all.plots$maxN$prev_10$BAY_LAG + ggtitle("E") + 
  plot_layout(design = hlay)
dev.off()
