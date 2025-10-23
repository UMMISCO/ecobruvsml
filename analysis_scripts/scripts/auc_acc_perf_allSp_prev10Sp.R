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

