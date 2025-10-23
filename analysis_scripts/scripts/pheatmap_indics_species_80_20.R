library(ggVennDiagram)
library(ggplot2)
library(readr) 
library(dplyr)
library(reshape2)

# load result list of 80_20 analysis (Indval/terinter/bininter)
load("/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/all_analysis_results_80_20.Rda")

names(alldatalist.video_80_20) <- ifelse(names(alldatalist.video_80_20)=="seed_10", "seed_10", gsub("seed_", "seed_0", names(alldatalist.video_80_20)))

# load indics species result table for analysis_80_20
indics_species.df <- read.csv(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/indics_species_analysis_results_80_20.csv")

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
indics_species_IN_OFF.df <- indics_species.df[indics_species.df$comparison=="INSHORE_OFFSHORE",]

# rename bin and ter with predomics_ter and predomics_bin
indics_species_IN_OFF.df[indics_species_IN_OFF.df == 'terinter'] <- 'predomics_ter'
indics_species_IN_OFF.df[indics_species_IN_OFF.df == 'bininter'] <- 'predomics_bin'

# get only species at revalence 10%

indics_species_IN_OFF_prev_10.df <- indics_species_IN_OFF.df[indics_species_IN_OFF.df$prev_rate=="prev_10",]

#Main loops to generate the Venn diagrams
for (i in unique(indics_species_IN_OFF_prev_10.df$seed)) {
  for (j in unique(indics_species_IN_OFF_prev_10.df$prev_rate)) {
    for (k in unique(indics_species_IN_OFF_prev_10.df$data)) {
      for (l in unique(indics_species_IN_OFF_prev_10.df$source)) {
        # Create Venn list for method
        venn_list_by_method[[i]][[j]][[k]][[l]] <- unique(indics_species_IN_OFF_prev_10.df$feature[
          (indics_species_IN_OFF_prev_10.df$seed == i) &
            (indics_species_IN_OFF_prev_10.df$prev_rate == j) &
            (indics_species_IN_OFF_prev_10.df$data == k) &
            (indics_species_IN_OFF_prev_10.df$source == l)
        ])
        # Create Venn list for method
        venn_list_by_source[[i]][[j]][[l]][[k]] <- venn_list_by_method[[i]][[j]][[k]][[l]]
      }
      # Generate Venn diagram by method
      venn_by_method[[i]][[j]][[k]] <- generate_venn_plot(venn_list_by_method[[i]][[j]][[k]], title_text = paste(i,j,k, sep="/"))
    }
  }
}

for (i in unique(indics_species_IN_OFF_prev_10.df$seed)) {
  for (j in unique(indics_species_IN_OFF_prev_10.df$prev_rate)) {
    for (k in unique(indics_species_IN_OFF_prev_10.df$source)) {
      # Generate Venn diagram by source
      venn_by_source[[i]][[j]][[k]] <- generate_venn_plot(venn_list_by_source[[i]][[j]][[k]], title_text = paste(i,j,k, sep="/"), set_size = 2.5)
    }
  }
}

input_data= list()

# get input data for each seed for prev_0 and prev_10 in INSHORE_OFFSHORE
for (i in unique(indics_species_IN_OFF_prev_10.df$seed)) {
  for (j in unique(indics_species_IN_OFF_prev_10.df$prev_rate)) {
    for (k in c("INSHORE_OFFSHORE")) {
      input_data[[i]][[j]][[k]]= alldatalist.video_80_20[[i]][[j]][["terinter"]][["predout.maxn"]][[k]][["Comp_data"]]
    }
  }
}

#names(input_data) <- gsub("seed_", "seed_0", names(input_data))

# merge indval, bininter and terinter species

#put all data frames into list
# indval_species= seed3_prev0_indics_species_IN_OFF_indval$feature
# bininter_species= seed3_prev0_indics_species_IN_OFF_pred_bin$feature
# terinter_species= seed3_prev0_indics_species_IN_OFF_pred_ter$feature
# seed3_prev0_maxN_list <- list(indval_species, bininter_species, terinter_species)      
# 
# seed3_prev0_bin <- list(seed3_prev0_indics_species_IN_OFF_indval, seed3_prev0_indics_species_IN_OFF_pred_bin, seed3_prev0_indics_species_IN_OFF_pred_ter)

#merge all data frames together
# seed3_prev0_merged_indval_bin_ter_presAbs <- Reduce(function(x, y) inner_join(x, y, by=c("feature","IsIndSp","comparison", "data", "seed", "prev_rate")), seed3_prev0_bin) 
# Reduce(intersect, list(a,b,c))

common_indics_species <- list()
for (i in unique(indics_species_IN_OFF_prev_10.df$seed)) {
  for (j in unique(indics_species_IN_OFF_prev_10.df$prev_rate)) {
    for (k in unique(indics_species_IN_OFF_prev_10.df$data)) {
      # Generate Venn diagram by source
      common_indics_species[[i]][[j]][[k]] <- Reduce(intersect, venn_list_by_method[[i]][[j]][[k]])
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

# pheatmap of indicator species by seed separated by train/test in Inshore-Offshore based on abundace
pheat_map.plots[["maxN"]]

# pheatmap of common indicator species by seed separated by train/test in Inshore-Offshore based on abundace
pheat_map_common.plots[["maxN"]]