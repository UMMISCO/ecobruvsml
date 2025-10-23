library(predomics)
library(dplyr)
library(ggalt)

load(file = "/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/all_analysis_results_80_20.Rda")

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
load(file = "/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/rf_analysis_results_80_20.Rda")

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

hlay="
A
A
B
C
C
D
"
p1 + 
  p1_dumbell + 
  p2 + 
  p2_dumbell + 
  plot_layout(design = hlay)



summary(auc_prev_10.dly.df$mean.value[!auc_prev_10.dly.df$method %in% c("BestModel_bininter","BestModel_terinter")])
summary(auc_prev_10.dly.df$mean.feats[!auc_prev_10.dly.df$method %in% c("BestModel_bininter","BestModel_terinter")])


summary(auc_prev_10.dly.df$mean.value[auc_prev_10.dly.df$method %in% c("BestModel_bininter","BestModel_terinter")])
summary(auc_prev_10.dly.df$mean.feats[auc_prev_10.dly.df$method %in% c("BestModel_bininter","BestModel_terinter")])



###########

# heatmaps train/test features

###########
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
          strip.text.y = element_text(angle = 360),
          axis.text.x = element_blank(),
          #axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

pdf(file="/data/projects/aime/analyses/bruvs/Figure5/Figure5.pdf", h=8, w=14)
hlay="
KKA
KKB
KKC
LLD
LLE
MMF
MMG
MMH
NNI
NNJ
"
x.plots$seed_01 + ggtitle("C") + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  p1 + ggtitle("A") + 
  p1_dumbell + theme(axis.text.y = element_blank()) +
  p2 + ggtitle("B") + 
  p2_dumbell + theme(axis.text.y = element_blank()) + 
  plot_layout(design = hlay)
dev.off()

pdf(file="/data/projects/aime/analyses/bruvs/Figure5/Figure5_hmapOnly.pdf", h=8, w=14)
hlay="
##A
##B
##C
##D
##E
##F
##G
##H
##I
##J
"
x.plots$seed_01 + ggtitle("C") + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  x.plots$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()) + 
  plot_layout(design = hlay)
dev.off()


# pheatmap of common indicator species by seed separated by train/test in Inshore-Offshore based on abundace
pheat_map_common.plots[["maxN"]]


hlay1="
A
A
B
"
g1 <- wrap_plots(p1,p1_dumbell+ theme(axis.text.y = element_blank()), design = hlay1)
g2 <- wrap_plots(p2,p2_dumbell+ theme(axis.text.y = element_blank()), design = hlay1)
g3 <-wrap_plots(x.plots$seed_01 + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                x.plots$seed_02 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                x.plots$seed_03 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_04 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_05 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_06 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_07 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_08 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_09 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()), 
                x.plots$seed_10 + theme(strip.text.x = element_blank(), strip.background.x = element_rect(fill = NA, linetype = 0), legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), strip.background.y = element_rect(fill = NA, linetype = 0), strip.text.y = element_blank()),
                ncol=1)
# save the updated version of Figure 5
pdf(file="/data/projects/aime/analyses/bruvs/Figure5_Update_EK/Figure5_Best.pdf", h=8, w=16)
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

