library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(ggvenn) # to ddisplay venn diagram with color set

dir.create("/data/projects/aime/analyses/bruvs/FigureS3/")

# load indval results
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/Indval_output_data/Indval_all_analyses_overall_data_prev_10.Rda")

# load bininter results
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/bininter_output_data/bininter_Predomics_all_analyses_overall_data_prev_10.Rda")

### get object for adonis results in all comarison from each algorithm
# Rename variables loaded to diffenriate it with terinter vars
adonis_pred_bininter.bin <- adonis_pred.bin; adonis_pred_bininter.maxn <- adonis_pred.maxn
predout_bininter.bin <- predout.bin; predout_bininter.bin.sub <- predout.bin.sub 
predout_bininter.maxn <- predout.maxn; predout_bininter.maxn.sub <- predout.maxn.sub
# remove loaded variables to avoid override
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

# load terinter results
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/terinter_output_data/terinter_Predomics_all_analyses_overall_data_prev_10.Rda")
# Rename variables loaded to diffenriate it with terinter vars
adonis_pred_terinter.bin <- adonis_pred.bin; adonis_pred_terinter.maxn <- adonis_pred.maxn
predout_terinter.bin <- predout.bin; predout_terinter.bin.sub <- predout.bin.sub
predout_terinter.maxn <- predout.maxn; predout_terinter.maxn.sub <- predout.maxn.sub
# remove loaded variables to avoid override
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

# initialize a list of adonis objects  
adonis.list <- list(ind.bin= adonis.bin, ind.maxn= adonis.maxn, pred_bin.bin= adonis_pred_bininter.bin, pred_bin.maxn= adonis_pred_bininter.maxn, pred_ter.bin= adonis_pred_terinter.bin, pred_ter.maxn= adonis_pred_terinter.maxn)

# Initialize an empty list to store the results
combined_adonis.list <- list()

# Loop through each element of the adonis.list
for (i in names(adonis.list)) {
  # Append each sub-list to the combined_adonis.list
  combined_adonis.list[[i]] <- do.call(rbind, adonis.list[[i]])
}

# Combine all the sub-lists into one dataframe by row
adonis.df <- do.call(rbind, combined_adonis.list)

# # Rename the "Pr..F." column of adonis.df
colnames(adonis.df)[colnames(adonis.df) == "Pr..F."] <- "pvalue"

# Replace all occurrences of the value 'allSpecies_terinter' or 'allSpecies_bininter' with 'allSpecies' in adonis.df
adonis.df$source[adonis.df$source %in% c('allSpecies_terinter', 'allSpecies_bininter')] <- 'allSpecies'

# computer the ratio between r2 of individual sites comparison and r2 of the community
library(dplyr)
adonis.df <- adonis.df %>%
  group_by(comparison, data) %>%
  dplyr::mutate(
    ratio = round(R2 / R2[source == "allSpecies"])
  ) %>%
  ungroup()

## get object for feature importance comparison of indics species from Indval/Predomics
# remove pval column from indvalout.bin.sub and indvalout.maxn.sub
indvalout.bin.sub <- indvalout.bin.sub[ , !(names(indvalout.bin.sub) %in% "pval")]
indvalout.maxn.sub <- indvalout.maxn.sub[ , !(names(indvalout.maxn.sub) %in% "pval")]

# Replace all occurrences of the value 'presAbs' with 'pres/abs' in indvalout.bin.sub
#indvalout.bin.sub[indvalout.bin.sub == 'presAbs'] <- 'pres/abs'

# initialize a list for all indicator species 
indics_species.list <- list(indval.bin =indvalout.bin.sub, indval.maxn= indvalout.maxn.sub,
                            bininter.bin= predout_bininter.bin.sub, bininter.maxn= predout_bininter.maxn.sub,
                            terinter.bin= predout_terinter.bin.sub, terinter.maxn= predout_terinter.maxn.sub)

# Combine all the sub-lists into one dataframe by row
indics_species.list <- do.call(rbind, indics_species.list)
# Remove row names of indics_species.list
rownames(indics_species.list) <- NULL

# Plot permanova results
adonis.df$colour <- ifelse(adonis.df$pvalue<=0.001,"blue","red")

# remove row of allspecies community before plotting
adonis.df_filt <- adonis.df[adonis.df$source != "allSpecies", ]

## Visualization of predomics results
library(predomics)
## redefine getImportanceFeaturesFBMobjects from predomics package
getImportanceFeaturesFBMobjects <- function (clf_res, X, y, verbose = TRUE, filter.cv.prev = 0, 
                                             scaled.importance = FALSE, k_penalty = 0, k_max = 0) 
{
  mode <- NULL
  if (!(predomics::isExperiment(clf_res))) {
    stop("analyzeLearningFeatures: please provide a valid experiment results!")
  }
  if (clf_res$classifier$params$objective == "cor") {
    mode <- "regression"
  }
  else {
    mode <- "classification"
  }
  if (!is.null(mode)) {
    cat(paste("... Estimating mode: ", mode, "\n"))
  }
  else {
    stop("analyzeImportanceFeatures: mode not founding stopping ...")
  }
  pop <- modelCollectionToPopulation(clf_res$classifier$models)
  if (verbose) 
    print(paste("There are", length(pop), "models in this population"))
  pop <- selectBestPopulation(pop = pop, p = 0.05, k_penalty = k_penalty, k_max = k_max)
  if (verbose) 
    print(paste("There are", length(pop), "models in this population after selection of the best"))
  # if (length(pop) == 1) {
  #   stop("analyzeImportanceFeatures: only one model after filtering. Plot can not be built... returing empty handed.")
  # }
  pop.df <- populationToDataFrame(pop = pop)
  pop.noz <- listOfModelsToDenseCoefMatrix(clf = clf_res$classifier, 
                                           X = X, y = y, list.models = pop)
  if (verbose) 
    print(paste("Pop noz object is created with", nrow(pop.noz), 
                "features and", ncol(pop.noz), "models"))
  fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf_res$classifier)
  pop.noz <- fa$pop.noz
  pop.noz <- data.frame(pop.noz)
  pop.noz$feature <- rownames(pop.noz)
  pop.noz <- melt(pop.noz)
  pop.noz$learner <- unlist(lapply(strsplit(as.character(pop.noz$variable), 
                                            split = "_"), function(x) {
                                              x[1]
                                            }))
  pop.noz$learner <- paste(pop.noz$learner, clf_res$classifier$params$language, 
                           sep = ".")
  pop.noz$model <- as.character(unlist(lapply(strsplit(as.character(pop.noz$variable), 
                                                       split = "_"), function(x) {
                                                         x[2]
                                                       })))
  pop.noz$value <- factor(pop.noz$value, levels = c(-1, 0, 
                                                    1))
  pop.noz$value <- droplevels(pop.noz$value)
  lr <- list(clf_res)
  names(lr) <- paste(clf_res$classifier$learner, clf_res$classifier$params$language, 
                     sep = ".")
  feat.import <- mergeMeltImportanceCV(list.results = lr, filter.cv.prev = filter.cv.prev, 
                                       min.kfold.nb = FALSE, learner.grep.pattern = "*", nb.top.features = nrow(X), 
                                       feature.selection = NULL, scaled.importance = scaled.importance, 
                                       make.plot = TRUE)
  if (is.null(feat.import)) {
    stop("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
  }
  if (mode == "regression") {
    featPrevPlot <- plotPrevalence(features = rownames(X), 
                                   X, y = NULL)
  }
  else {
    featPrevPlot <- plotPrevalence(features = rownames(X), 
                                   X, y)
  }
  effSizes.df <- computeEffectSizes(X = X, y = y, mode = mode)
  outlist <- list(featprevFBM = pop.noz, featImp = feat.import$summary, 
                  effectSizes = effSizes.df, featPrevGroups = featPrevPlot$data)
  class(outlist) <- "listFBMfeatures"
  return(outlist)
}

alldatacomparison <- list()

for (i in unique(adonis.df$comparison)){
  alldatacomparison[[i]][['terinter']] = list(maxN= predout_terinter.maxn[[i]], presAbs=predout_terinter.bin[[i]])
  alldatacomparison[[i]][['bininter']] = list(maxN= predout_bininter.maxn[[i]], presAbs=predout_bininter.bin[[i]])
}

video.fmbObj =list()
fbm.plots =list()
library(reshape2)
for (i in unique(adonis.df$comparison)) {
  for (j in c("maxN", "presAbs")){
    for (k in c("bininter", "terinter")) {
      pred.obj=alldatacomparison[[i]][[k]][[j]]
      res_clf=pred.obj$fit
      X= pred.obj$Comp_data$X
      y= pred.obj$Comp_data$y
      video.fmbObj[[i]][[j]][[k]] <- getImportanceFeaturesFBMobjects(clf_res = res_clf, X = t(X), y = y, verbose = TRUE, k_penalty = 0, k_max=0, filter.cv.prev =0)
    }
    fbm.plots[[i]][[j]] <- plotImportanceFeaturesFBMobjects(FBMobjList = video.fmbObj[[i]][[j]], verbose = TRUE, nb.top.features = 100, makeplot = FALSE)
  }
}

# plots FBM + FI + CliffDelta + FeatPrev
#wrap_plots(wrap_elements(fbm.plots[['maxN']]), wrap_elements(fbm.plots[['bin']]), nrow = 2) + plot_annotation(tag_levels = "A")
fbm.plot= wrap_plots(wrap_elements(fbm.plots$BAR_BAY$maxN),
                     wrap_elements(fbm.plots$BAR_BAY$presAbs),
                     wrap_elements(fbm.plots$BAR_LAG$maxN),
                     wrap_elements(fbm.plots$BAR_LAG$presAbs),
                     wrap_elements(fbm.plots$BAY_LAG$maxN),
                     wrap_elements(fbm.plots$BAY_LAG$presAbs), nrow = 3, ncol=2) + plot_annotation(tag_levels = "A")

# save FBM + FI + CliffDelta + FeatPrev for pairwise classification to a pdf
pdf(file="/data/projects/aime/analyses/bruvs/FigureS3/FigS3.pdf", h=15, w=25)
fbm.plot
dev.off()
