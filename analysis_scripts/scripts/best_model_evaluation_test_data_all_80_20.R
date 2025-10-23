library(predomics)
library(dplyr)

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
  #geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "AUC") +
  facet_grid(.~source) +
  # ylim(c(0.6,1.05))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = after_stat(p.signif)), label.y = c(1.01, 1.02, 1.03,1.04)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12), strip.text.x = element_text(size=12))+
  # scale_size(range = c(0.3,1)) +
  scale_y_continuous(breaks = seq(0.3, 1, 0.1))

# select only perf for prev_10
acc_prev_10.df <- acc.df[acc.df$prev_rate=='prev_10', ]

# plot acc perf for prev_10
acc_prev_10.plot <- ggplot(data=acc_prev_10.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  #geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "ACC") +
  facet_grid(.~source) +
  # ylim(c(0.6,1.05))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = after_stat(p.signif)), label.y = c(1.01, 1.02, 1.03,1.04)) +
  theme(legend.text = element_text(size = 12),
        # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8),
        strip.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  # scale_size(range = c(0.3,1)) +
  scale_y_continuous(breaks = seq(0.3, 1, 0.1))

# select only perf for prev_10
f1_prev_10.df <- f1.df[f1.df$prev_rate=='prev_10', ]

# plot f1 perf for prev_10
f1_prev_10.plot <- ggplot(data=f1_prev_10.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  #geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "f1") +
  facet_grid(.~source) +
  # ylim(c(0.6,1.05))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = after_stat(p.signif)), label.y = c(1.01, 1.02, 1.03,1.04)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), 
        # legend.text = element_text(size = 12),
        legend.position = "none",
        # strip.text.x = element_text(size=12),
        strip.text = element_blank(),)+
  # scale_size(range = c(0.3,1)) +
  scale_y_continuous(breaks = seq(0.3, 1, 0.1))

# combine ACC and f1 for prev_10 in one plot
wrap_plots(acc_prev_10.plot, f1_prev_10.plot, nrow = 2)

# plot all metrics performances of RF
# my_comparisons <- list(c("Indval", "predomics_bin"), c("Indval", "predomics_ter"), c("predomics_bin", "predomics_ter"))

ggplot(data=auc.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  #geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "AUC") +
  facet_grid(prev_rate~source) +
  ylim(c(0.5,1.25))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = after_stat(p.signif)), label.y = c(1.05, 1.1, 1.15, 1.2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))#+
  #scale_size(range = c(1,5), limits = c(1,80))

ggplot(data=acc.df, aes(x= comparison, y= value, color=method, label=featNum))+
  geom_boxplot(alpha = 0.6)+
  #geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "ACC") +
  facet_grid(prev_rate~source) +
  ylim(c(0.3,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))#+
  #scale_size(range = c(1,5), limits = c(1,80))

ggplot(data=f1.df, aes(x= comparison, y= value, color=method))+
  geom_boxplot(alpha = 0.6)+
  #geom_point(aes(size=featNum), position = position_jitterdodge(dodge.width = 0.7))+
  labs(x = "Comparison", y = "f1") +
  facet_grid(prev_rate~source) +
  ylim(c(0.3,1.5))+
  theme_bw() +
  stat_compare_means(aes(group = method, label = paste0("p = ", after_stat(p.format))), label.y = c(1.1, 1.2, 1.3, 1.4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 8), legend.text = element_text(size = 12))#+
  #scale_size(range = c(1,5), limits = c(1,80))

# plot number of species by method on inshore-offshore for prev_10
auc_prev_10_IN_OFF.df= auc_prev_10.df[auc_prev_10.df$comparison=="INSHORE_OFFSHORE",]

# barplot of feature number by model and source
ggplot(data=auc_prev_10_IN_OFF.df, aes(x=seed, y=featNum,fill=method)) +
  geom_bar(stat="identity", position="dodge")+
  facet_grid(.~source)+
  geom_text(aes(label=featNum), vjust=-1.1, size=3, position=position_dodge(1.2))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# lineplot of feature number by model and source
ggplot(auc_prev_10_IN_OFF.df, aes(x = seed, y = featNum, color = method, group = method)) +
  geom_line(size = 1) +  # Add lines for trends
  geom_point(size = 2) + # Add points for clarity
  facet_wrap(~source) + # Separate panels for maxN and presAbs
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
    text = element_text(size = 12), # Larger text for clarity
    strip.text = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "bottom" # Place legend below the plot
  )

# library(dplyr)
# 
# all_indics_species.df$source= as.factor(all_indics_species.df$source)
# 
# list1.indval <- all_indics_species.df %>%
#   filter(prev_rate == "prev_0", data == "maxN", IsIndSp==1, source=="indval")%>%
#   pull(feature)
# 
# 
# list2.indval <- all_indics_species.df %>%
#   filter(prev_rate == "prev_0", data == "pres/abs", IsIndSp==1, source=="indval")%>%
#   pull(feature)
# 
# list3.indval <- all_indics_species.df %>%
#   filter(prev_rate == "prev_10", data == "maxN", IsIndSp==1, source=="indval")%>%
#   pull(feature)
# 
# list4.indval <- all_indics_species.df %>%
#   filter(prev_rate == "prev_10", data == "pres/abs", IsIndSp==1, source=="indval")%>%
#   pull(feature)
# 
# list1.predomics_bin <- all_indics_species.df %>%
#   filter(prev_rate == "prev_0", data == "maxN", IsIndSp==1, source=="predomics_bin")%>%
#   pull(feature)
# 
# 
# list2.predomics_bin <- all_indics_species.df %>%
#   filter(prev_rate == "prev_0", data == "pres/abs", IsIndSp==1, source=="predomics_bin")%>%
#   pull(feature)
# 
# list3.predomics_bin <- all_indics_species.df %>%
#   filter(prev_rate == "prev_10", data == "maxN", IsIndSp==1, source=="predomics_bin")%>%
#   pull(feature)
# 
# list4.predomics_bin <- all_indics_species.df %>%
#   filter(prev_rate == "prev_10", data == "pres/abs", IsIndSp==1, source=="predomics_bin")%>%
#   pull(feature)
# 
# all_indics_species.df$source= as.factor(all_indics_species.df$source)
# 
# list1.predomics_ter <- all_indics_species.df %>%
#   filter(prev_rate == "prev_0", data == "maxN", IsIndSp==1, source=="predomics_ter")%>%
#   pull(feature)
# 
# 
# list2.predomics_ter <- all_indics_species.df %>%
#   filter(prev_rate == "prev_0", data == "pres/abs", IsIndSp==1, source=="predomics_ter")%>%
#   pull(feature)
# 
# list3.predomics_ter <- all_indics_species.df %>%
#   filter(prev_rate == "prev_10", data == "maxN", IsIndSp==1, source=="predomics_ter")%>%
#   pull(feature)
# 
# list4.predomics_ter <- all_indics_species.df %>%
#   filter(prev_rate == "prev_10", data == "pres/abs", IsIndSp==1, source=="predomics_ter")%>%
#   pull(feature)
# 
# indics.species.plot.list= list(prev0_maxN.indval= list1.indval, prev10_maxN.indval= list2.indval, prev0_presAbs.indval=list3.indval, prev10_presAbs.indval=list4.indval,
#                                prev0_maxN.predomics_bin= list1.predomics_bin, prev10_maxN.predomics_bin= list2.predomics_bin, prev0_presAbs.predomics_bin=list3.predomics_bin, prev10_presAbs.predomics_bin=list4.predomics_bin,
#                                prev0_maxN.predomics_ter= list1.predomics_ter, prev10_maxN.predomics_ter= list2.predomics_ter, prev0_presAbs.predomics_ter=list3.predomics_ter, prev10_presAbs.predomics_ter=list4.predomics_ter                               )
# 
# UpSetR::upset(fromList(indics.species.plot.list), nsets=12,
#       sets = c("prev0_maxN.indval", "prev10_maxN.indval", "prev0_presAbs.indval","prev10_presAbs.indval",
#                "prev0_maxN.predomics_bin", "prev10_maxN.predomics_bin", "prev0_presAbs.predomics_bin","prev10_presAbs.predomics_bin",
#                "prev0_maxN.predomics_ter", "prev10_maxN.predomics_ter", "prev0_presAbs.predomics_ter","prev10_presAbs.predomics_ter"),
#       keep.order = TRUE, order.by ="freq",
# 
# )
# 
# upset(indics.species.plot.list, nsets = 12, sets = names(indics.species.plot.list))
# 
# grep("\\.indval$", strings, value = TRUE)
# 
# library("UpSetR")
