library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(readr)

# load indval results
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/Indval_output_data/Indval_all_analyses_overall_data_prev_0.Rda")

# load bininter results
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/bininter_output_data/bininter_Predomics_all_analyses_overall_data_prev_0.Rda")

### get object for adonis results in all comarison from each algorithm
# Rename variables loaded to diffenriate it with terinter vars
adonis_pred_bininter.bin <- adonis_pred.bin; adonis_pred_bininter.maxn <- adonis_pred.maxn
predout_bininter.bin <- predout.bin; predout_bininter.bin.sub <- predout.bin.sub 
predout_bininter.maxn <- predout.maxn; predout_bininter.maxn.sub <- predout.maxn.sub
# remove loaded variables to avoid override
rm(adonis_pred.bin, adonis_pred.maxn, predout.bin, predout.bin.sub, predout.maxn, predout.maxn.sub)

# load terinter results
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/terinter_output_data/terinter_Predomics_all_analyses_overall_data_prev_0.Rda")
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

# # Remove specified columns: "Df", "SumOfSqs", "F"
# adonis.df <- adonis.df[, !colnames(adonis.df) %in% c("Df", "SumOfSqs", "F")]
# # Rename the "Pr..F." column of adonis.df
colnames(adonis.df)[colnames(adonis.df) == "Pr..F."] <- "pvalue"
# Replace all occurrences of the value 'pres/abs' with 'bin' in adonis.df
# adonis.df[adonis.df == 'pres/abs'] <- 'bin'
# Replace all occurrences of the value 'allSpecies_terinter' or 'allSpecies_bininter' with 'allSpecies' in adonis.df
adonis.df$source[adonis.df$source %in% c('allSpecies_terinter', 'allSpecies_bininter')] <- 'allSpecies'
# # Remove rows where "source" is equal to any value in the list c("allSpecies_terinter", "allSpecies_bininter")
# adonis.df <- adonis.df[!adonis.df$source %in% c("allSpecies_terinter", "allSpecies_bininter"), ]

# # remove duplicated rows of adonis.df
# adonis.df <- adonis.df %>%
#   distinct(comparison, data, source, .keep_all = TRUE)

# computer the ratio between r2 of individual sites comparison and r2 of the community
library(dplyr)
adonis.df <- adonis.df %>%
  group_by(comparison, data) %>%
  mutate(
    ratio = round(R2 / R2[source == "allSpecies"], digits = 3)
  ) %>%
  ungroup()


## get object for feature importance comparison of indics species from Indval/Predomics
# remove pval column from indvalout.bin.sub and indvalout.maxn.sub
indvalout.bin.sub <- indvalout.bin.sub[ , !(names(indvalout.bin.sub) %in% "pval")]
indvalout.maxn.sub <- indvalout.maxn.sub[ , !(names(indvalout.maxn.sub) %in% "pval")]

# Reorder columns of indvalout.bin.sub and indvalout.maxn.sub
# indvalout.bin.sub <- indvalout.bin.sub[, c("feature","featureImportance","class","IsIndSp","prevalence","comparison","source","data")]
# indvalout.maxn.sub <- indvalout.maxn.sub[, c("feature","featureImportance","class","IsIndSp","prevalence","comparison","source","data")]


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

# get permanova results for only INSHORE8OFFSHORE comparison
adonis.df_IN_OFF = adonis.df[adonis.df$comparison=="INSHORE_OFFSHORE",]

# remove row of allspecies community before plotting
adonis.df_filt <- adonis.df[adonis.df$source != "allSpecies", ]

# filter only INSHORE_OFFSHORE comparison
adonis.df_filt_IN_OFF = adonis.df_filt[adonis.df_filt$comparison=="INSHORE_OFFSHORE",]

adonis.df_filt_IN_OFF$colour <- ifelse(adonis.df_filt_IN_OFF$pvalue<=0.001,"blue","red")

# barplot of pemonova by r2
adonis.df_filt_IN_OFF$source= as.factor(adonis.df_filt_IN_OFF$source)
#pdf(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-D/permanova_all_results_updated_final.pdf", h=6, w=7)
ggplot(data=adonis.df_IN_OFF, aes(x=data, y=R2, fill=source, color=colour)) +
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(label=ifelse(pvalue <= 0.001, "**" , ifelse(pvalue <= 0.05, "*" , " " ))), vjust=-0.3, size=5,position=position_dodge(.9))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c('blue', 'red'), labels =  c( 'pval <= 0.001', 'pval > 0.001 & pval <=0.05')) +
  labs(color = "Significativity") 
#dev.off()

# barplot of pemonova for r2 of each comparison under r2 of INSHORE_OFFSHORE
#pdf(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-D/permanova_ratio_all_results_final.pdf", h=6, w=7)
# ggplot(data=adonis.df_filt, aes(x=comparison, y=ratio, color=colour)) +
#   geom_bar(stat="identity", fill="lightblue")+
#   geom_text(aes(label=ratio), vjust=-0.7, size=3.5,nudge_y=-0.0003)+ 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   labs(color = "Significativity") +
#   scale_color_manual(values = c('blue', 'red'), labels =  c( 'pval <= 0.001', 'pval > 0.001 & pval <=0.05')) +
#   scale_y_continuous(name="r2_comp/r2_community", limits=c(0, 6)) +
#   facet_grid(data~source)
#dev.off()
# plot R2, ratio and Significativity by comparison
ggplot(data=adonis.df, aes(x=comparison, y=R2,fill=source)) +
     geom_bar(stat="identity", position = "dodge")+
     geom_text(aes(label=ifelse(pvalue <= 0.001, "**" , ifelse(pvalue <= 0.05, "*" , " " ))), vjust=-0.03, size=5, position=position_dodge(.9))+
     geom_text(aes(label=ratio), vjust=2.5, size=3.5, position=position_dodge(.9))+
     theme(axis.text.x = element_text(angle = 45, hjust = 1))+
     facet_grid(data~.)

# plot ratio R2_indics_species/R2_community on inshore/offshore
ggplot(data=adonis.df_IN_OFF, aes(x=data, y=ratio, fill=source, color=colour)) +
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(label=ifelse(pvalue <= 0.001, "**" , ifelse(pvalue <= 0.05, "*" , " " ))), vjust=-0.1, size=5, position=position_dodge(.9))+
  geom_text(aes(label=ratio), color = "black", vjust=2.5, size=3.5, position=position_dodge(.9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(color = "Significativity") +
  scale_color_manual(values = c('blue', 'red'), labels =  c( 'pval <= 0.001', 'pval > 0.001 & pval <=0.05')) +
  scale_y_continuous(name="r2_comp/r2_community", limits=c(0, 2))


# plot number of indics species by method on inshore-offshore
indSpNumber.plot <- ggplot(data=adonis.df_filt_IN_OFF, aes(x=source, y=features, fill=source)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_manual(values = c("#648C67","#8195B2","#CC3837")) +
  ylab("Number of indicator species") +
  xlab("Methods") +
  facet_grid(.~data) +
  theme_bw()+
  geom_text(aes(label=features), vjust=-0.9, size=5, position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 12),
        legend.position = "none") +
  scale_y_continuous(limits=c(0, 47)) +
  scale_x_discrete(labels=c("Indval", "Pred_Terinter", "Pred_Bininter"))

# filter indics species list by inshore_offshore
all_species_IN_OFF.df <- indics_species.list[indics_species.list$comparison=="INSHORE_OFFSHORE",]

# filter only indicator species
all_indics_species_IN_OFF.df <- all_species_IN_OFF.df[all_species_IN_OFF.df$IsIndSp==1,]

# rename bin and ter with predomics_ter and predomics_bin
all_species_IN_OFF.df[all_species_IN_OFF.df == 'terinter'] <- 'Pred_Terinter'
all_species_IN_OFF.df[all_species_IN_OFF.df == 'bininter'] <- 'Pred_Bininter'
all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df == 'terinter'] <- 'Pred_Terinter'
all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df == 'bininter'] <- 'Pred_Bininter'

# Required library
library(ggVennDiagram)
library(ggvenn)
library(patchwork)

# Initialize lists
venn_list_by_method <- list()
venn_by_method <- list()
venn_list_by_source <- list()
venn_by_source <- list()

# Helper function for generating Venn diagrams
generate_venn_plot <- function(venn_data, title_text, set_size = 3) {
  #ggVennDiagram(venn_data, label= "count", set_size = set_size, label_size = 4) +
  ggvenn(venn_data, fill_color = c("#648C67","#CC3837","#8195B2"), stroke_size = 0.5, set_name_size = 4,
         show_percentage = FALSE, stroke_alpha = 1)+
    ggtitle(title_text) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(color = "black", size = 9, margin = margin(0, 10, 0, 0), hjust = 0.5),
      legend.position = "none"
    ) + guides(fill = guide_legend(title = "IndicSpecies")) +
    #scale_fill_distiller(palette = "RdBu", limits=c(0,15)) +
    scale_x_continuous(expand = expansion(mult = 0.25))
  
}

# library(ggalluvial)

# # Plot using ggalluvial
# ggplot(all_species_IN_OFF.df, aes(x = source, stratum = IsIndSp, alluvium = data, y = feature)) +
#   scale_x_discrete(expand = c(.1, .1)) +
#   geom_alluvium() +
#   geom_stratum(alpha = .5) +
#   geom_text(stat = "stratum", size = 3) +
#   theme(legend.position = "none") +
#   ggtitle("vaccination survey responses at three points in time")
# 
# # Corrected alluvial plot
# alluvium_data= to_alluvia_form(all_species_IN_OFF.df, key=source, value =IsIndSp, id=data)
# ggplot(all_species_IN_OFF.df, aes(x = source, stratum = IsIndSp, alluvium = data, y = feature)) +
#   geom_alluvium(aes(fill = IsIndSp), width = 0.2) +
#   geom_stratum(alpha = 0.5, fill = "gray") +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
#   scale_x_discrete(expand = c(0.1, 0.1)) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   ggtitle("Vaccination survey responses at three points in time")

# Main loops to generate the Venn diagrams
for (i in unique(all_indics_species_IN_OFF.df$data)) {
  for (j in unique(all_indics_species_IN_OFF.df$source)) {
    # Create Venn list for method
    venn_list_by_method[[i]][[j]] <- unique(all_indics_species_IN_OFF.df$feature[
      (all_indics_species_IN_OFF.df$data == i) &
        (all_indics_species_IN_OFF.df$source == j)
    ])
    # Create Venn list for source
    venn_list_by_source[[j]][[i]] <- venn_list_by_method[[i]][[j]]
  }
  
  # Generate Venn diagram by method
  ifelse(i=='maxN',
         venn_by_method[[i]] <- generate_venn_plot(venn_list_by_method[[i]], title_text = "Abundance (maxN)"),
         venn_by_method[[i]] <- generate_venn_plot(venn_list_by_method[[i]], title_text = "Presence/absence (pres/abs)"))
}

for (k in unique(all_indics_species_IN_OFF.df$source)) {
  # Generate Venn diagram by method
  venn_by_source[[k]] <- generate_venn_plot(venn_list_by_source[[k]], title_text = k, set_size = 2.5)
}
# plot common indicator species in different method across abundance and pres/abs
# venn_by_source_merged <- (venn_by_source[['indval']])/(venn_by_source[['pred_bin']]+venn_by_source[['pred_ter']]) + plot_layout(heights = c(1, 1)) 

# plot common indicator species in abundance and presence/absence across method
# venn_by_method_merged <- venn_by_method[['maxN']] + venn_by_method[['pres/abs']]

# get common species names for maxN table
common.sp.maxn <- intersect(intersect(venn_list_by_method[['maxN']]$indval, venn_list_by_method[['maxN']]$Pred_Bininter), venn_list_by_method[['maxN']]$Pred_Terinter)

# get common species names for pres/abs table
common.sp.bin <- intersect(intersect(venn_list_by_method[['pres/abs']]$indval, venn_list_by_method[['pres/abs']]$Pred_Bininter), venn_list_by_method[['pres/abs']]$Pred_Terinter)

# Define the plots
plot_maxN_method <- venn_by_method[['maxN']]
plot_pres_abs_method <- venn_by_method[['pres/abs']]

# Combine plots with a blank spacer between them
combined_venn_by_method <- plot_maxN_method + plot_pres_abs_method +
  plot_layout(widths = c(1.2, 1.2)) # Adjust the 0.1 to make the spacer narrower or wider

# plot number of indics species with Venn diagram
indSpNumber.plot/(combined_venn_by_method) + plot_annotation(tag_levels = list(c("A", "B")))

## Visualization of predomics Models performances (Terinter, bininter)

###########

## First integrated visualization= performance of models in training/testing steps

##############

# combine all results in one list
alldatalist.video.IN_OFF = list()
alldatalist.video.IN_OFF$bininter = list(maxN= predout_bininter.maxn[["INSHORE_OFFSHORE"]], bin=predout_bininter.bin[["INSHORE_OFFSHORE"]])
alldatalist.video.IN_OFF$terinter = list(maxN= predout_terinter.maxn[["INSHORE_OFFSHORE"]], bin=predout_terinter.bin[["INSHORE_OFFSHORE"]])

# Build a list to save the needed datasets for plotting
library(reshape2)
library(ggpubr)
digest.list <- list()
for (m in names(alldatalist.video.IN_OFF)) {
  for (i in c("maxN", "bin")) {
    iobj <- alldatalist.video.IN_OFF[[m]][[i]][["digest"]]
    #iterate over acc and auc tables in cv$scores (source of boxplot visualization in digest plots)
    for(t in c("generalization.auc","generalization.acc"))
    {
      #print(t)
      #Get the data
      tdf <- iobj$cv$scores[[t]]
      #exclude the na
      tdf <- tdf[rowSums(!is.na(tdf))>0,]
      #melt + process for plotting
      tdf <- data.frame(parsimony = rownames(tdf), tdf)
      tdf <- melt(tdf, id.vars = "parsimony")
      #Add the comparison (the i)
      tdf$comparison <- "INSHORE_OFFSHORE"
      #Add the variable (here habitat)
      tdf$source <- i
      # Add the algorithm
      tdf$algorithm <- m
      #Save
      digest.list[[t]][[i]][[m]] <- tdf
    }
    #iterate over the empirical performance data of best model in each prediction task (source of lineplot visualization in digest plots)
    for(t in c("accuracy_","auc_","recall_","precision_"))
    {
      #print(t)
      tdf <- data.frame(value=iobj$best$scores[[t]], 
                        parsimony=names(iobj$best$scores[[t]]),
                        comparison="INSHORE_OFFSHORE",
                        source=i,
                        algorithm= m)
      #Save
      digest.list[[t]][[i]][[m]] <- tdf
    }
  }
}

#For each sub-element in digest.list, do the rbind of the dataframes of each binary prediction task
for (n in names(digest.list)) {
  digest.list[[n]] <- lapply(digest.list[[n]], function(x){do.call("rbind", x)})
}
#For each element in digest.list, do the rbind of the dataframes of each binary prediction task
digest.list <- lapply(digest.list, function(x){do.call("rbind", x)})

#Build the visualization plots, first boxplots of acc/auc across training/test CV data
digest.plots <- list()
for(k in c("generalization.auc","generalization.acc"))
{
  kdf <- digest.list[[k]]
  #Factorize the parsimony variable
  kdf$algorithm <- factor(kdf$algorithm)
  # kdf$comparison <- factor(kdf$comparison,
  #                         levels = levels(kdf$parsimony)[order(as.numeric(gsub("k_","", levels(kdf$parsimony))))])
  #Do the integrated plot (here the facet_grid(.~comparison) command stratifies the plot by comparison)
  # plot by source
  digest.plots[[k]] <- ggplot(data = kdf, aes(y = value, x = algorithm)) +
    geom_point(aes(color = "dark"), position = position_jitterdodge(dodge.width = 0.9), size = 1, alpha = 0.5)+
    geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width = 0.9), alpha = 0.5) +
    ylab(gsub("^.*\\.","", k)) +
    xlab("predomics_method") +
    ggtitle(gsub("\\..*$","", k)) +
    ylim(c(0,1.25)) +
    facet_grid(.~source) +
    stat_compare_means(method = "wilcox.test",
                       label.y = 1.15)+
    theme_bw() +
    # geom_hline(yintercept = unique(maj.class),col = "gray", linetype = "dashed") +
    theme(legend.position = "bottom",
          legend.direction = "horizontal") +
    guides(colour = "none") #+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6)) # +
  # theme(strip.text.x=element_text(angle=90, size=8))
}
#Second, lineplots of empirical performances
for(i in c("accuracy_","auc_"))
{
  idf <- digest.list[[i]]
  #Factorize the parsimony variable
  idf$parsimony <- factor(idf$parsimony)
  idf$parsimony <- factor(idf$parsimony,
                          levels = levels(idf$parsimony)[order(as.numeric(gsub("k_","", levels(idf$parsimony))))])
  ibest <- list()
  #get the data for best model in each case as dataframe
  for(c in unique(idf$algorithm))
  {
    for (j in c("maxN", "bin")) {
      best.model.obj <- alldatalist.video.IN_OFF[[c]][[j]]$digest$best$model
      ibest[[j]][[c]] <- data.frame(learner=best.model.obj$learner,
                               language=best.model.obj$language,
                               variable=best.model.obj[[i]],
                               variable.source=i,
                               k=length(best.model.obj$indices_),
                               algorithm=c, 
                               source=j)
    }
  }
  
  #For each element in digest.list, do the rbind of the dataframes of each binary prediction task
  ibest <- lapply(ibest, function(x){do.call("rbind", x)})
  
  #merge individual dataframes in single dataframe (rbind passed to do.call function + list ibest)
  ibest <- do.call("rbind",ibest)
  #Add the labels that we want to visualize corresponding to the best model of each prediction task
  ibest$label <- paste("L:", ibest$learner,"_", ibest$language,"; F:", signif(ibest$variable,2), "; K=", ibest$k, sep = "")
  ibest$x <- (length(unique(idf$parsimony))+1)/2
  ibest$y <- seq(1,nrow(ibest),1)/20
  #do the integrated plot (colour lines and points by comparisons + add the infos for the best model in each comparison stored in the ibest table) + save
  digest.plots[[i]] <- ggplot(data = idf, aes(x = parsimony, y = value, group = 1)) +
    geom_line(aes(group=algorithm,color = algorithm)) +
    geom_point(size = 2, alpha = 1) +
    geom_text(data=ibest, aes(x=x,y=y, label=label, colour=algorithm), size=3, inherit.aes = FALSE) +
    ylab(i) +
    xlab("Model parsimony") +
    ggtitle("Emprical performance") +
    ylim(c(0,1)) +
    theme_bw() +
    facet_grid(.~source) +
    # geom_hline(yintercept = mean(maj.class), col = "gray", linetype = "dashed") +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.title = element_blank())+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))
}

#Integrated visualization as in the output of digest function but now integrating all pairwise comparisons across levels of habitat variable
# pdf(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/models_performances_pred_video_data_pa.pdf", h=4, w=6)
grid.arrange(digest.plots$generalization.auc,
  digest.plots$generalization.acc,
  digest.plots$auc_,
  digest.plots$accuracy_,
  ncol=2, nrow=2)

# dev.off()

# save auc and acc performances in generalization and empirical
auc_acc_all_species.list <- list()
auc_acc_all_species.list$generalization.auc <- digest.list$generalization.auc
auc_acc_all_species.list$empirical.auc <- digest.list$auc_
auc_acc_all_species.list$generalization.acc <- digest.list$generalization.acc
auc_acc_all_species.list$empirical.acc <- digest.list$accuracy_

# match dataframe structure
for (i in names(auc_acc_all_species.list)) {
  if( i %in% c("generalization.auc", "generalization.acc")){
    auc_acc_all_species.list[[i]] = auc_acc_all_species.list[[i]][, !names(auc_acc_all_species.list[[i]]) %in% "variable"]
  }
  auc_acc_all_species.list[[i]][["metric"]] <- i
  # reorder the columns
  auc_acc_all_species.list[[i]] <- auc_acc_all_species.list[[i]][, c("metric", "value","parsimony","comparison","source","algorithm")]
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

# rbind auc and acc perf data of auc_acc_all_species.list

auc_acc_all_species.df <- flatten_and_bind(auc_acc_all_species.list)
rownames(auc_acc_all_species.df) <- NULL

# add a column for dataset
auc_acc_all_species.df$dataset <- "all_species_community"

# save model performances
# write.csv(auc_acc_all_species.df, file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/auc_acc_all_species.csv")

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

video.fmbObj =list()
fbm.plots =list()
for (i in c("maxN", "bin")){
  for (j in names(alldatalist.video.IN_OFF)) {
    pred.obj=alldatalist.video.IN_OFF[[j]][[i]]
    res_clf=pred.obj$fit
    X= pred.obj$Comp_data$X
    y= pred.obj$Comp_data$y
    video.fmbObj[[i]][[j]] <- getImportanceFeaturesFBMobjects(clf_res = res_clf, X = t(X), y = y, verbose = TRUE, k_penalty = 0, k_max=0, filter.cv.prev =0)
  }
  fbm.plots[[i]] <- plotImportanceFeaturesFBMobjects(FBMobjList = video.fmbObj[[i]], verbose = TRUE, nb.top.features = 100, makeplot = FALSE)
}
# plots FBM + FI + CliffDelta + FeatPrev
wrap_plots(wrap_elements(fbm.plots[['maxN']]), wrap_elements(fbm.plots[['bin']]), nrow = 2) + plot_annotation(tag_levels = "A")

# plot prevalence distribution

all_indics_species_IN_OFF.df$source <- as.factor(all_indics_species_IN_OFF.df$source)
# distribution of density indics species vs others
ggplot(all_indics_species_IN_OFF.df, aes(x = prevalence, color=source, group=source)) +
  geom_density(alpha = 0.6) +
  labs(x = "Prevalence", y = "Density") +
  scale_x_continuous(limits=c(0, 100))+
  #scale_color_manual(values = c("pred_bin" = "#FF6666", "indval" = "#87CEEB", "pred_ter" = "#90EE90"))+
  facet_grid(.~data)+
  theme_bw() 


# plot comparison between Indval/terinter/bininter by prevalence
comparison= combn(x = as.character(unique(all_species_IN_OFF.df$IsIndSp)), m = 2, simplify = FALSE)
library(ggpubr)
# factorize IsIndSp
all_species_IN_OFF.df$IsIndSp= as.factor(all_species_IN_OFF.df$IsIndSp)

# prevalence distribution betweeen species that are indicator or not
prev_distribution <- ggplot(all_species_IN_OFF.df, aes(x=source, y=prevalence, color=IsIndSp))+
  geom_boxplot(alpha=0.6)+
  labs(x = "Source", y = "Prevalence") +
  facet_grid(. ~data) +
  theme_bw()+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.y = c(61, 55, 52))

# scatter plot prevalence vs feature importance
# factorize IsIndSp
all_indics_species_IN_OFF.df$IsIndSp= as.factor(all_indics_species_IN_OFF.df$IsIndSp)
all_indics_species_IN_OFF_indval.df= all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df$source=="indval",]
# Rename the feat import in all_species_indval.df to indicator value column of adonis.df
colnames(all_indics_species_IN_OFF_indval.df)[colnames(all_indics_species_IN_OFF_indval.df) == "featureImportance"] <- "indicator_value"

all_indics_species_IN_OFF_bininter.df= all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df$source=="Pred_Bininter",]

all_indics_species_IN_OFF_terinter.df= all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df$source=="Pred_Terinter",]

## scatter plot feat import vs prevalence
library(ggpubr)

ggplot(all_indics_species_IN_OFF_indval.df, aes(x=indicator_value, y=prevalence))+
  geom_point()+
  stat_cor(method = "spearman")+
  labs(x = "indicator_value", y = "Prevalence") +
  facet_grid(.~data) + ggtitle("Indval")+
  theme_bw()

ggplot(all_indics_species_IN_OFF_bininter.df, aes(x=featureImportance, y=prevalence))+
  geom_point()+
  stat_cor(method = "spearman")+
  labs(x = "featureImportance", y = "Prevalence") +
  facet_grid(.~data) + ggtitle("Bininter")+
  theme_bw()

ggplot(all_indics_species_IN_OFF_terinter.df, aes(x=featureImportance, y=prevalence))+
  geom_point()+
  stat_cor(method = "spearman")+
  labs(x = "featureImportance", y = "Prevalence") +
  facet_grid(.~data) + ggtitle("Terinter")+
  theme_bw()

all_species_IN_OFF.df$IsIndSp= as.factor(all_species_IN_OFF.df$IsIndSp)
all_species_IN_OFF_indval.df= all_species_IN_OFF.df[all_species_IN_OFF.df$source=="indval",]
# Rename the feat import in all_species_indval.df to indicator value column of adonis.df
colnames(all_species_IN_OFF_indval.df)[colnames(all_species_IN_OFF_indval.df) == "featureImportance"] <- "indicator_value"

all_species_IN_OFF_bininter.df= all_species_IN_OFF.df[all_species_IN_OFF.df$source=="Pred_Bininter",]

all_species_IN_OFF_terinter.df= all_species_IN_OFF.df[all_species_IN_OFF.df$source=="Pred_Terinter",]
## indval vs bininter
merged_indval_bin = merge(all_species_IN_OFF_indval.df, all_species_IN_OFF_bininter.df, by=c("feature","comparison","data"), all.x = TRUE)

## indval vs bininter vs terinter
merged_indval_bin_ter = merge(merged_indval_bin, all_species_IN_OFF_terinter.df, by=c("feature","comparison","data"), all.x = TRUE)
#merged_indval_ter_clean = na.omit(merged_indval_ter)
# View(inval_terinter_species)

indval_species = merged_indval_bin_ter[merged_indval_bin_ter$IsIndSp.x==1,]
bininter_species = merged_indval_bin_ter[merged_indval_bin_ter$IsIndSp.y==1,]
terinter_species = merged_indval_bin_ter[merged_indval_bin_ter$IsIndSp==1,]

# get indics_species from indval and bininter
indval_bininter_species= rbind(indval_species, bininter_species)
indval_bininter_species$color= ifelse(indval_bininter_species$IsIndSp.x==1 & indval_bininter_species$IsIndSp.y==0,
                                      -1, 
                                      ifelse(indval_bininter_species$IsIndSp.x==1 & indval_bininter_species$IsIndSp.y==1, 0, 1))
indval_bininter_species$color <- factor(indval_bininter_species$color, levels=c(-1,0,1), labels = c("indVal-only","indVal-Pred","Pred-only"))

# get indics_species from indval and terinter
indval_terinter_species= rbind(indval_species, terinter_species)
indval_terinter_species$color= ifelse(indval_terinter_species$IsIndSp.x==1 & indval_terinter_species$IsIndSp==0,
                                      -1, 
                                      ifelse(indval_terinter_species$IsIndSp.x==1 & indval_terinter_species$IsIndSp==1, 0, 1))
indval_terinter_species$color <- factor(indval_terinter_species$color, levels=c(-1,0,1), labels = c("indVal-only","indVal-Pred","Pred-only"))

# correlation between indval and bininter
indval.feat1 = ggplot(unique(indval_bininter_species), aes(x=featureImportance.x, y=indicator_value)) + 
  geom_point(aes(color=color)) + 
  stat_cor(method = "spearman", aes(label = after_stat(p.label)), label.y = 1) +
  labs(x= "feature_importance", y= "indicator_value", title = "Indval vs bininter", color = "species_source")+
  facet_grid(. ~ data) +
  theme_bw()+
  theme(legend.position = "none")

indval.feat1

# correlation between indval and terinter
indval.feat2 = ggplot(unique(indval_terinter_species), aes(x=featureImportance.y, y=indicator_value)) + 
  geom_point(aes(color=color)) +
  stat_cor(method = "spearman", aes(label = after_stat(p.label)), label.y = 1) +
  labs(x= "feature_importance", y= "indicator_value", title = "Indval vs terinter", color = "species_source")+
  facet_grid(. ~ data) +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

indval.feat2

## combine 4 figures in one for article: number of indics species by method, venn diagram of indics species shared by all methods,
##  prevalence distribution between indics and non-indics species, correlation between indicator_value and feature_importance

indSpNumber.plot_wrap   <- patchwork::wrap_elements(indSpNumber.plot)
prev_distribution_wrap   <- patchwork::wrap_elements(prev_distribution)
indval.feat_import <- patchwork::wrap_elements(indval.feat1 + indval.feat2)
venn_by_source     <- patchwork::wrap_elements(venn_by_method[['maxN']] + venn_by_method[['pres/abs']])

# combine plots of number of indics species, Venn diagram, prevalence distribution and correlation indicator_value vs feat importance
result.plot <- patchwork::wrap_plots(indSpNumber.plot_wrap,
                                     venn_by_source,
                                     prev_distribution_wrap,
                                     indval.feat_import) + plot_annotation(tag_levels = "A")

# get indicator species for all methods
indval_bin_ter_species.df= merged_indval_bin_ter
# rename IsIndSp based on the method
colnames(indval_bin_ter_species.df)[colnames(indval_bin_ter_species.df) == "IsIndSp.x"] <- "IsIndSp.indval"
colnames(indval_bin_ter_species.df)[colnames(indval_bin_ter_species.df) == "IsIndSp.y"] <- "IsIndSp.pred_bin"
colnames(indval_bin_ter_species.df)[colnames(indval_bin_ter_species.df) == "IsIndSp"] <- "IsIndSp.pred_ter"

# Annotate common or specific species by method
indval_bin_ter_species.df$IndSp_by_method <- ifelse(indval_bin_ter_species.df$IsIndSp.indval == 0 & indval_bin_ter_species.df$IsIndSp.pred_bin == 0 & indval_bin_ter_species.df$IsIndSp.pred_ter == 1, 1,
                                                           ifelse(indval_bin_ter_species.df$IsIndSp.indval == 0 & indval_bin_ter_species.df$IsIndSp.pred_bin == 1 & indval_bin_ter_species.df$IsIndSp.pred_ter == 0, 2, 
                                                                  ifelse(indval_bin_ter_species.df$IsIndSp.indval == 1 & indval_bin_ter_species.df$IsIndSp.pred_bin == 0 & indval_bin_ter_species.df$IsIndSp.pred_ter == 0, 3,
                                                                         ifelse(indval_bin_ter_species.df$IsIndSp.indval == 1 & indval_bin_ter_species.df$IsIndSp.pred_bin == 0 & indval_bin_ter_species.df$IsIndSp.pred_ter == 1, 4, 
                                                                                ifelse(indval_bin_ter_species.df$IsIndSp.indval == 1 & indval_bin_ter_species.df$IsIndSp.pred_bin == 1 & indval_bin_ter_species.df$IsIndSp.pred_ter == 0, 5, 
                                                                                       ifelse(indval_bin_ter_species.df$IsIndSp.indval == 1 & indval_bin_ter_species.df$IsIndSp.pred_bin == 1 & indval_bin_ter_species.df$IsIndSp.pred_ter == 1, 6, NA)
                                                                                )
                                                                         )
                                                                  )
                                                           )
)

# Remove NAs
indval_bin_ter_species.df <- indval_bin_ter_species.df[!is.na(indval_bin_ter_species.df$IndSp_by_method), ]

# add labels to common species levels
indval_bin_ter_species.df$IndSp_by_method <- factor(indval_bin_ter_species.df$IndSp_by_method, levels=c(1,2,3,4,5,6), labels = c("Pred_ter-only","Pred_bin-only","indval-only","indval-Pred_ter","indval-Pred_bin","indval-Pred_bin_ter"))

# plot the indicator species with their prevalence by labels
ggplot(indval_bin_ter_species.df, aes(x=prevalence, y=feature, fill=IndSp_by_method)) +
     geom_bar(stat="identity") +
     guides(fill=guide_legend(reverse=TRUE))+
     facet_grid(data~.)+
     theme_bw()

ggplot(indval_bin_ter_species.df, aes(x=prevalence, y=feature)) +
  geom_segment(aes(yend=feature), xend=0, colour="grey50") +
  geom_point(size=3, aes(colour=data)) +
  facet_grid(IndSp_by_method~., scales="free", space="free") +
  #scale_colour_brewer(palette="Set1", limits=c("NL","AL"), guide="none")+
  theme_bw() +
  theme(axis.text.x = element_blank())

## plot using alluvial diagram
library(ggalluvial)
ggplot(data = indval_bin_ter_species.df,
       aes(axis1 = feature,   # First variable on the X-axis
           axis2 = IndSp_by_method, 
           y = prevalence)) +
  geom_alluvium(aes(fill = prevalence)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("feature", "IndSp_by_method"),
                   expand = c(0.005, 0.15)) +
  theme_void()


## NMDS (Non-metric multidimensional distance scaling) on community species, indval species, terinter species and bininter species on Bray-curtis distance
library(vegan)
library(dplyr)

# load the data
load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/video_data_object.rda")

# compute nmds on community data
community.maxn.X= predout_bininter.maxn$INSHORE_OFFSHORE$Comp_data$X
community.bin.X= predout_bininter.bin$INSHORE_OFFSHORE$Comp_data$X
nmds.community.maxn= metaMDS(community.maxn.X, distance = "bray")
nmds.community.bin= metaMDS(community.bin.X, distance = "jaccard")

## select data for indval species and compute nmds
indval_data.maxn = filter(indval_species, data == "maxN")
indval_data.bin = filter(indval_species, data == "pres/abs")
indval_data.maxn.X = community.maxn.X[, colnames(community.maxn.X) %in% indval_data.maxn$feature]
indval_data.bin.X = community.bin.X[, colnames(community.bin.X) %in% indval_data.bin$feature]
nmds.indval.maxn= metaMDS(indval_data.maxn.X, distance = "bray")
nmds.indval.bin= metaMDS(indval_data.bin.X, distance = "jaccard")
#nmds.indval.scores = as.data.frame(vegan::scores(nmds.indval, "species"))
    
# compute nmds on bininter species data
bininter_data.maxn = filter(bininter_species, data == "maxN")
bininter_data.bin = filter(bininter_species, data == "pres/abs")
bininter_data.maxn.X = community.maxn.X[, colnames(community.maxn.X) %in% bininter_data.maxn$feature]
bininter_data.bin.X = community.bin.X[, colnames(community.bin.X) %in% bininter_data.bin$feature]
# remove sample with no species
bininter_data.maxn.X= bininter_data.maxn.X[rowSums(bininter_data.maxn.X>0)>0,]
bininter_data.bin.X= bininter_data.bin.X[rowSums(bininter_data.bin.X>0)>0,]
# compute nmds
nmds.bininter.maxn= metaMDS(bininter_data.maxn.X, distance = "bray")
nmds.bininter.bin= metaMDS(bininter_data.bin.X, distance = "jaccard")
#nmds.indval.scores = as.data.frame(vegan::scores(nmds.indval, "species"))

# compute nmds on terinter species data
terinter_data.maxn = filter(terinter_species, data == "maxN")
terinter_data.bin = filter(terinter_species, data == "pres/abs")
terinter_data.maxn.X = community.maxn.X[, colnames(community.maxn.X) %in% terinter_data.maxn$feature]
terinter_data.bin.X = community.bin.X[, colnames(community.bin.X) %in% terinter_data.bin$feature]
# remove sample with no species
terinter_data.maxn.X= terinter_data.maxn.X[rowSums(terinter_data.maxn.X>0)>0,]
terinter_data.bin.X= terinter_data.bin.X[rowSums(terinter_data.bin.X>0)>0,]
# compute nmds
nmds.terinter.maxn= metaMDS(terinter_data.maxn.X, distance = "bray")
nmds.terinter.bin= metaMDS(terinter_data.bin.X, distance = "jaccard")

# store NMDS data in a list
nmds_data =list()
nmds_data[["maxN"]][["community"]] = nmds.community.maxn
nmds_data[["maxN"]][["indval"]] = nmds.indval.maxn
nmds_data[["maxN"]][["terinter"]] = nmds.terinter.maxn
nmds_data[["maxN"]][["bininter"]] = nmds.bininter.maxn
nmds_data[["presAbs"]][["community"]] = nmds.community.bin
nmds_data[["presAbs"]][["indval"]] = nmds.indval.bin
nmds_data[["presAbs"]][["terinter"]] = nmds.terinter.bin
nmds_data[["presAbs"]][["bininter"]] = nmds.bininter.bin


# plotting NMDS of each set of species

# function for ellipsess 
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
library(ggplot2)
library(ggpubr)
library(ggrepel)

nmds.plots = list()
for (i in names(nmds_data)){
  for (j in c("community","indval","bininter", "terinter")) {
    nmds = nmds_data[[i]][[j]]
    nmds.sites = data.frame(nmds$points)
    nmds.sites$Habitat = sm$sample_info$Habitat[rownames(sm$sample_info) %in% rownames(nmds.sites)]
    nmds.sites$Station <- rownames(nmds.sites)
    nmds.sites$Site <- factor(sm$sample_info$Site[match(nmds.sites$Station, sm$sample_info$Station)], levels = c("ABAR","MBAR","ALAG","MLAG","ABAY","MBAY"))
    nmds.sites$Site <- as.character(nmds.sites$Site)
    
    nmds.sites$Site[nmds.sites$Site %in% "ABAR"] <- "Barrier"
    nmds.sites$Site[nmds.sites$Site %in% "MBAR"] <- "Barrier"
    nmds.sites$Site[nmds.sites$Site %in% "ALAG"] <- "Lagoon"
    nmds.sites$Site[nmds.sites$Site %in% "MLAG"] <- "Lagoon"
    nmds.sites$Site[nmds.sites$Site %in% "ABAY"] <- "Bay"
    nmds.sites$Site[nmds.sites$Site %in% "MBAY"] <- "Bay"
    # nmds.sites$Site[nmds.sites$Site %in% "ABAR"] <- "Abor?? Barrier"
    # nmds.sites$Site[nmds.sites$Site %in% "MBAR"] <- "M'b??r?? Barrier"
    # nmds.sites$Site[nmds.sites$Site %in% "ALAG"] <- "Abor?? Lagoon"
    # nmds.sites$Site[nmds.sites$Site %in% "MLAG"] <- "M'b??r?? Lagoon"
    # nmds.sites$Site[nmds.sites$Site %in% "ABAY"] <- "Abor?? Bay"
    # nmds.sites$Site[nmds.sites$Site %in% "MBAY"] <- "M'b??r?? Bay"
    #nmds.sites$Site <- factor(nmds.sites$Site, levels = c("Abor?? Barrier", "M'b??r?? Barrier","Abor?? Lagoon", "M'b??r?? Lagoon","Abor?? Bay","M'b??r?? Bay" ))
    unique(nmds.sites$Site)
    
    nmds.species = data.frame(nmds$species)
    nmds.species$Species = rownames(nmds.species)
    
    dfEllSite <- data.frame() #sets up a data frame before running the function.
    # Loop over each level of `Site` in `nmds.sites`
    for (g in levels(nmds.sites$Site)) {
      
      # Subset the data for the current site
      site_data <- nmds.sites[nmds.sites$Site == g, ]
      
      # Check if there is data available for this site
      if (nrow(site_data) > 1) {
        
        # Remove any rows with NA values in MDS1 or MDS2 to avoid errors
        site_data <- site_data[complete.cases(site_data$MDS1, site_data$MDS2), ]
        
        # Calculate the covariance and center if we have non-zero data points
        if (nrow(site_data) > 1) {
          cov_ellipse <- veganCovEllipse(
            cov.wt(cbind(site_data$MDS1, site_data$MDS2), wt = rep(1 / nrow(site_data), nrow(site_data)))$cov,
            center = c(mean(site_data$MDS1), mean(site_data$MDS2))
          )
          
          # Combine the ellipse data with the Site identifier
          dfEllSite <- rbind(dfEllSite, cbind(as.data.frame(cov_ellipse), Site = g))
        }
      }
    }
    
    siteScoresMeanPts <- aggregate(nmds.sites[ ,c("MDS1", "MDS2")],
                                   list(group = nmds.sites$Site), mean)
    
    nmds.plots[[i]][[j]] <- ggplot(data = nmds.sites, aes(x=MDS1, y=MDS2, colour=Site)) +
      geom_point(shape=10, size=2) + # add the point markers
      geom_path(data = dfEllSite, aes(x = V1, y = V2, group = Site), linewidth = 1.5) +
      #geom_text(data=spScores,aes(x=MDS1,y=MDS2,label= species),alpha=0.5, color = "grey50", size = 3.5) +  # add the species labels
      ggtitle(paste(i,j, sep = "|"))+
      geom_text_repel(data=nmds.species,
                      aes(x=MDS1,y=MDS2,label=Species, # hjust=0.5*(1-sign(MDS1)),vjust=0.5*(1-sign(MDS2))
                      ),
                      color="grey50", segment.color = "grey80", min.segment.length	= 0.5, alpha=0.5, size=2.5, fontface = "italic",
                      max.overlaps = 15) + #, box.padding = 0.5, force = 0.1, force_pull = 2
      
      #scale_colour_manual(values=c("Abor?? Barrier" = "#0072B2", "Abor?? Bay" = "#D55E00", "Abor?? Lagoon" = "#009E73", "M'b??r?? Barrier" = "#56B4E9", "M'b??r?? Bay" = "#E69F00", "M'b??r?? Lagoon" = "#F0E442")) +
      scale_colour_manual(values=c("Abor?? Barrier" = "#0072B2", "Abor?? Bay" = "#D55E00", "Abor?? Lagoon" = "#009E73", "M'b??r?? Barrier" = "#56B4E9", "M'b??r?? Bay" = "#E69F00", "M'b??r?? Lagoon" = "#F0E442")) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", colour = "black"),
            legend.justification = c(0,1),
            legend.position= "bottom",
            legend.direction = "horizontal",
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            legend.background = element_blank()
      )
  }
  
}
library(patchwork)
library(gridExtra)
# # combine all nmds plots
nmds_abund <- grid.arrange(nmds.plots$maxN$community,
                           nmds.plots$maxN$indval, 
                           nmds.plots$maxN$terinter, 
                           nmds.plots$maxN$bininter,
                           ncol=2, nrow=2, widths = c(2.7,2.7))
nmds_presAbs <- grid.arrange(nmds.plots$presAbs$community,
                           nmds.plots$presAbs$indval, 
                           nmds.plots$presAbs$terinter, 
                           nmds.plots$presAbs$bininter,
                           ncol=2, nrow=2, widths = c(2.7,2.7))

nmds.combined = wrap_plots(nmds_abund,nmds_presAbs) + plot_annotation(tag_levels = 'A')


##############################################################################################################################################
# filter indics species list by inshore_offshore
adonis.list1 <- adonis.df_filt_IN_OFF
adonis.list1$prev_rate= "prev_0"

adonis.list2 <- adonis.df_filt_IN_OFF
adonis.list2$prev_rate= "prev_10"

adonis.df= bind_rows(list(adonis.list1, adonis.list2))


# filter indics species list by inshore_offshore
all_indics_species.list1 <- indics_species.list[indics_species.list$comparison=="INSHORE_OFFSHORE",]
all_indics_species.list1$prev_rate= "prev_0"

# filter only indicator species
indics_species.list1 <- all_indics_species.list1[all_indics_species.list1$IsIndSp==1,]


# filter indics species list by inshore_offshore
all_indics_species.list2 <- indics_species.list[indics_species.list$comparison=="INSHORE_OFFSHORE",]
all_indics_species.list2$prev_rate= "prev_10"

# filter only indicator species
indics_species.list2 <- all_indics_species.list2[all_indics_species.list2$IsIndSp==1,]

# merge dataframes
all_indics_species.df= bind_rows(list(all_indics_species.list1,all_indics_species.list2))

indics_species.df= bind_rows(list(indics_species.list1, indics_species.list2))

# plot comparison between Indval/terinter/bininter by prevalence
comparison= combn(x = as.character(unique(all_indics_species.df$IsIndSp)), m = 2, simplify = FALSE)
library(ggpubr)
# factorize IsIndSp
all_indics_species.df$IsIndSp= as.factor(all_indics_species.df$IsIndSp)

ggplot(all_indics_species.df, aes(x=source, y=prevalence, color=IsIndSp))+
  geom_boxplot(alpha=0.6)+
  labs(x = "Source", y = "Prevalence") +
  facet_grid(prev_rate~data) +
  theme_bw()+
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.y = c(69, 67, 65))+
theme(axis.text.x = element_text(angle = 45, hjust = 1))

# rename 
all_indics_species.df[all_indics_species.df == 'terinter'] <- 'predomics_ter'
all_indics_species.df[all_indics_species.df == 'bininter'] <- 'predomics_bin'

# plot number of indics species by method on inshore-offshore
indSpNumber.plot <- ggplot(data=adonis.df, aes(x=data, y=features,fill=source)) +
  geom_bar(stat="identity", position = "dodge")+
  geom_text(aes(label=features), vjust=-1.5, size=5, position=position_dodge(.9))+
  facet_grid(prev_rate~.) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits=c(0, 50))

# filter indics species list by inshore_offshore
all_species_IN_OFF.df <- indics_species.list[indics_species.list$comparison=="INSHORE_OFFSHORE",]

# filter only indicator species
all_indics_species_IN_OFF.df <- all_species_IN_OFF.df[all_species_IN_OFF.df$IsIndSp==1,]



# rename bin and ter with predomics_ter and predomics_bin
all_species_IN_OFF.df[all_species_IN_OFF.df == 'terinter'] <- 'predomics_ter'
all_species_IN_OFF.df[all_species_IN_OFF.df == 'bininter'] <- 'predomics_bin'
all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df == 'terinter'] <- 'pred_ter'
all_indics_species_IN_OFF.df[all_indics_species_IN_OFF.df == 'bininter'] <- 'pred_bin'

# Required library
library(ggVennDiagram)
library(patchwork)

# Initialize lists
venn_list_by_method <- list()
venn_by_method <- list()
venn_list_by_source <- list()
venn_by_source <- list()

# Helper function for generating Venn diagrams
generate_venn_plot <- function(venn_data, title_text, set_size = 3) {
  ggVennDiagram(venn_data, set_size = set_size, label_size = 2) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(color = "black", size = 9, margin = margin(10, 0, 10, 0)),
      legend.position = "bottom"
    ) + guides(fill = guide_legend(title = "IndicSpecies")) +
    scale_fill_distiller(palette = "RdBu", limits=c(0,40))
}

# Main loops to generate the Venn diagrams
for (i in unique(all_indics_species_IN_OFF.df$data)) {
  for (j in unique(all_indics_species_IN_OFF.df$source)) {
    # Create Venn list for method
    venn_list_by_method[[i]][[j]] <- unique(all_indics_species_IN_OFF.df$feature[
      (all_indics_species_IN_OFF.df$data == i) &
        (all_indics_species_IN_OFF.df$source == j)
    ])
    # Create Venn list for source
    venn_list_by_source[[j]][[i]] <- venn_list_by_method[[i]][[j]]
  }
  
  # Generate Venn diagram by method
  venn_by_method[[i]] <- generate_venn_plot(venn_list_by_method[[i]], title_text = i)
}

for (k in unique(all_indics_species_IN_OFF.df$source)) {
  # Generate Venn diagram by method
  venn_by_source[[k]] <- generate_venn_plot(venn_list_by_source[[k]], title_text = k, set_size = 2.5)
}
# plot common indicator species in different method across abundance and pres/abs
venn_by_source <- (venn_by_source[['indval']])/(venn_by_source[['pred_bin']]+venn_by_source[['pred_ter']]) + plot_layout(heights = c(1, 1)) 

# plot common indicator species in abundance and presence/absence across method
# venn_by_method <- venn_by_method[['maxN']] + venn_by_method[['pres/abs']]

# Define the plots
plot_maxN_method <- venn_by_method[['maxN']]
plot_pres_abs_method <- venn_by_method[['pres/abs']]

# Combine plots with a blank spacer between them
combined_venn_by_method <- plot_maxN_method + plot_pres_abs_method +
  plot_layout(widths = c(1.2, 1.2)) # Adjust the 0.1 to make the spacer narrower or wider

# plot number of indics species with Venn diagram
indSpNumber.plot/(combined_venn_by_method) + plot_annotation(tag_levels = list(c("A", "B")))


# import dplyr 
library(dplyr)
library(ggplot2)
library(ggpubr)

# script to evaluate auc and acc performances of predomics models in INSHORE/OFFSHORE for empirical and generalization for 5 cross-validation

# import dataframe of auc and acc performances for all species
auc_acc_all_species.df <- read.csv(file = "/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/auc_acc_all_species.csv")

# import dataframe of auc and acc performances for species with at least 10% of prevalence
auc_acc_species_prev_10.df <- read.csv(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/auc_acc_species_prev_10.csv")

# bind rows to get one dataframe
auc_acc.df <- dplyr::bind_rows(auc_acc_all_species.df, auc_acc_species_prev_10.df)

# filter only AUC performances
auc.df <- auc_acc.df[auc_acc.df$metric %in% c("generalization.auc", "empirical.auc"), ]

## Ensure source is a factor, ordered alphabetically, and add custom labels
source_labels <- c("maxN" = "maxN", "bin" = "pres/abs")  # Replace with actual values

auc.df$source <- factor(auc.df$source, 
                        levels = unique(auc.df$source), 
                        labels = source_labels[unique(auc.df$source)])

# Order the dataframe by the 'source' column in ascending order
auc_ordered.df <- auc.df[order(auc.df$source), ]

## violin plot of AUC performances of Predomics models on all_species_community vs species_filtered
auc_violin.plot <- ggplot(auc_ordered.df, aes(x = metric, y = value, fill = dataset)) +
  geom_violin(trim = TRUE, alpha = 0.5, position = position_dodge(width = 0.9), scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9), alpha = 0.7) +
  xlab("Metrics") +
  ylab("AUC value") +
  facet_grid(algorithm ~ source) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(angle = 0),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.7))

# plot AUC of performances of Predomics models on all_species_community vs species_filtered
auc_acc.plot <- ggplot(data = auc.df, aes(y = value, x = metric, fill=dataset)) +
  #geom_point(aes(color=value), position = position_jitterdodge(dodge.width = 0.2), size = 1, alpha = 0.5)+
  geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width = 0.9), alpha = 0.7) +
  xlab("Dataset") +
  ylab("metric_value") +
  #ylim(c(0,1.25)) +
  facet_grid(algorithm~source) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE)+
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        strip.text.x = element_text(angle = 0),
        axis.text.x = element_text(margin = margin(1,0.5,1,0.5,"cm"), angle = 30, vjust = 0.5, hjust=0.7))


