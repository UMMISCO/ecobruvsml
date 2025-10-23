library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(readr)
library(ggvenn)
library(ggVennDiagram)
library(ggpubr)
library(grid)
library(reshape2)
library(viridis)

dir.create("/data/projects/aime/analyses/bruvs/Figure3/")

objlist <- list()
# load indval results
x <- load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/Indval_output_data/Indval_all_analyses_overall_data_prev_10.Rda")
for(i in x)
{
  objlist[["indval"]][[i]] <- get(i)
}
rm(list = x)

# load bininter results
x <- load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/bininter_output_data/bininter_Predomics_all_analyses_overall_data_prev_10.Rda")
for(i in x)
{
  objlist[["bininter"]][[i]] <- get(i)
}
rm(list = x)

# load terinter results
x <- load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/terinter_output_data/terinter_Predomics_all_analyses_overall_data_prev_10.Rda")
for(i in x)
{
  objlist[["terinter"]][[i]] <- get(i)
}
rm(list = x)

mergedf.maxN <- rbind(objlist$indval$indvalout.maxn.sub[,-match("pval", colnames(objlist$indval$indvalout.maxn.sub))],
                      objlist$bininter$predout.maxn.sub,
                      objlist$terinter$predout.maxn.sub)

mergedf.bin <- rbind(objlist$indval$indvalout.bin.sub[,-match("pval", colnames(objlist$indval$indvalout.bin.sub))],
                     objlist$bininter$predout.bin.sub,
                     objlist$terinter$predout.bin.sub)

mergedf.all <- rbind(mergedf.bin, mergedf.maxN)
mergedf.all.bplot <- mergedf.all[mergedf.all$IsIndSp == 1,]
mergedf.all.bplot <- data.frame(table(mergedf.all.bplot$source, mergedf.all.bplot$data, mergedf.all.bplot$comparison))

colours.method <- c("#648C67","#8195B2","#CC3837")

plots.list <- list()
for(i in unique(mergedf.all.bplot$Var3))
{
  print(i)
  idf <- mergedf.all.bplot[mergedf.all.bplot$Var3 %in% i,]
  plots.list[[i]][["barplot"]] <- ggplot(idf, aes(x=Var2, y=Freq, fill=Var1)) + 
    geom_bar(stat="identity", position="dodge") + 
    geom_text(aes(label=Freq), vjust=0, position=position_dodge(.9)) + 
    scale_fill_manual(values = colours.method) + 
    ylab("Indicator Species") + 
    xlab("Data Source") + 
    labs(fill="Method") + 
    theme_bw()
}

for(i in unique(mergedf.all$comparison))
{
  print(i)
  idf <- mergedf.all[mergedf.all$comparison %in% i,]
  idf <- idf[idf$IsIndSp %in% 1,]
  i.listVenn <- list()
  for(j in unique(idf$data))
  {
    i.listVenn[[j]] <- list(bininter=idf$feature[idf$source %in% "bininter" & idf$data %in% j],
                            indval=idf$feature[idf$source %in% "indval" & idf$data %in% j],
                            terinter=idf$feature[idf$source %in% "terinter" & idf$data %in% j])
  }
  plots.list[[i]][["vennDiagrams_byMethod"]] <- lapply(i.listVenn, function(x){ggvenn(x,fill_color = c("#648C67","#8195B2","#CC3837"), show_stats = "c", stroke_size = 0.5, set_name_size = 3.5, stroke_alpha = 1, padding = 0.2) + 
      theme(plot.title = element_text(hjust = 0.5, vjust=0.5),
            # plot.margin = margin(0, 0.5,0 ,0 , "cm"),
            # theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 1, unit = "cm")),plot.margin = margin(2, 2, 2, 2, "mm"),
            # legend.title = element_text(color = "black", size = 9, margin = margin(0, 10, 0, 0), hjust = 0.5),
            legend.position = "none")  +
      scale_x_continuous(expand = expansion(mult = 0.01))})
  
  i.listVenn2 <- list()
  for(j in unique(idf$source))
  {
    i.listVenn2[[j]] <- list(maxN=idf$feature[idf$source %in% j & idf$data %in% "maxN"],
                             presabs=idf$feature[idf$source %in% j & idf$data %in% "pres/abs"])
  }
  plots.list[[i]][["vennDiagrams_byData"]] <- lapply(i.listVenn2, function(x){ggvenn(x, stroke_size = 0.5, set_name_size = 3.5, show_stats = "c", stroke_alpha = 1, padding = 0.2) + 
      ggtitle(x)+
      theme(plot.title = element_text(hjust = 0.5, vjust=1),
            # plot.margin = margin(0, 0.5,0 ,0 , "cm"),
            # legend.title = element_text(color = "black", size = 9, margin = margin(0, 10, 0, 0), hjust = 0.5),
            legend.position = "none")  +
      scale_x_continuous(expand = expansion(mult = 0.1))})
}

grid.arrange(plots.list$INSHORE_OFFSHORE$barplot,
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$maxN,
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$`pres/abs`, ncol=3)

adonis.bin.indval <- do.call("rbind", objlist$indval$adonis.bin)
adonis.bin.bininter <- do.call("rbind", objlist$bininter$adonis_pred.bin)
adonis.bin.terinter <- do.call("rbind", objlist$terinter$adonis_pred.bin)
adonis.bin <- rbind(adonis.bin.indval, adonis.bin.bininter, adonis.bin.terinter)

adonis.maxn.indval <- do.call("rbind", objlist$indval$adonis.maxn)
adonis.maxn.bininter <- do.call("rbind", objlist$bininter$adonis_pred.maxn)
adonis.maxn.terinter <- do.call("rbind", objlist$terinter$adonis_pred.maxn)
adonis.maxn <- rbind(adonis.maxn.indval, adonis.maxn.bininter, adonis.maxn.terinter)

adonis.alldf <- rbind(adonis.bin, adonis.maxn)
adonis.alldf$pval.cat <- ifelse(adonis.alldf$Pr..F.<0.05, "*","")
adonis.alldf <- adonis.alldf[!adonis.alldf$source %in% c("allSpecies_bininter","allSpecies_terinter"),]
adonis.alldf$source <- gsub("Species"," Species", adonis.alldf$source)
adonis.alldf$source <- factor(adonis.alldf$source, levels = c("all Species","terinter Species", "indval Species", "bininter Species"))


library(viridis)
ggplot(adonis.alldf, aes(x=source, y=comparison)) + 
  geom_point(aes(size=features, colour=R2)) + 
  scale_color_viridis() + 
  facet_grid(.~data) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

adonis.plots <- list()
for(i in unique(adonis.alldf$comparison))
{
  print(i)
  adonis.plots[[i]] <- ggplot(adonis.alldf[adonis.alldf$comparison %in% i,], aes(x=source, y=R2)) + 
    geom_bar(stat = "identity", aes(fill=features)) + 
    geom_text(aes(label=pval.cat)) + 
    scale_fill_viridis() + 
    facet_grid(data~.) + 
    coord_flip() +
    theme_bw()
}

do.call("grid.arrange",c(adonis.plots, ncol=2, nrow=2))

prevplots <- list()
for(i in unique(mergedf.all$comparison))
{
  print(i)
  idf <- mergedf.all[mergedf.all$comparison %in% i,]
  idf$IsIndSp <- factor(idf$IsIndSp, levels=c(0,1), labels=c("No","Yes"))
  prevplots[[i]] <- ggplot(idf, aes(x=source, y=prevalence)) + 
    geom_boxplot(aes(fill=IsIndSp), outlier.shape = NA, alpha=.5) + 
    geom_point(aes(colour=IsIndSp), position = position_dodge(width = .8)) + 
    facet_grid(.~data) + 
    scale_color_manual(values = c("black", "grey")) +
    scale_fill_manual(values = c("black", "grey")) +
    stat_compare_means(aes(group=IsIndSp), method = "wilcox", label = "p.signif") + 
    ylab("Species prevalence") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

corr.list <- list()
# boxplots.featImp <- list()
# for(i in unique(mergedf.all$comparison))
# {
#   print(i)
#   idf <- mergedf.all[mergedf.all$comparison %in% i,]
#   ##indval-bininter
#   idf1 <- idf[idf$source %in% c("indval","bininter"),]
#   idf1.dcast <- dcast(data = idf1, formula = feature+data~source, value.var="featureImportance")
#   idf1.dcast2 <- dcast(data = idf1, formula = feature+data~source, value.var="IsIndSp")
#   idf1.dcast2$class <- ifelse((idf1.dcast2$bininter==1 & idf1.dcast2$indval==1),"common",
#                               ifelse((idf1.dcast2$bininter==1 & idf1.dcast2$indval==0),"bininterOnly",
#                                      ifelse(idf1.dcast2$bininter==0 & idf1.dcast2$indval==1,"indvalOnly",NA)))
#   idf1.dcast2 <- idf1.dcast2[!is.na(idf1.dcast2$class),]
#   ilist.bininter <- list()
#   for(j in unique(idf1.dcast2$class))
#   {
#     print(j)
#     for(m in unique(idf1.dcast2$data))
#     {
#       print(m)
#       ##Get the features 
#       jfeats <- idf1.dcast2$feature[idf1.dcast2$class %in% j & idf1.dcast2$data %in% m]
#       ##Correlate feature importance metrics if >1
#       if(length(jfeats)>1)
#       {
#         jfeats.corr <- cor.test(idf1.dcast$bininter[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m], 
#                                 idf1.dcast$indval[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m], 
#                                 method="spearman")
#         jfeats.corr.df <- data.frame(comparison="indval-bininter",
#                                      feature.class=j,
#                                      dataSource=m,
#                                      rho=jfeats.corr$estimate,
#                                      pval=jfeats.corr$p.value,
#                                      features=length(jfeats))
#         ilist.bininter[[paste(j,m)]] <- jfeats.corr.df
#       }
#     }
#   }
#   ilist.bininter <- do.call("rbind", ilist.bininter)
#   ##boxplots featureImportance by feature class
#   idf1.dcast$mvar <- paste(idf1.dcast$feature, idf1.dcast$data)
#   colnames(idf1.dcast)[3:4] <- paste(colnames(idf1.dcast)[3:4],".FI", sep = "")
#   idf1.dcast2$mvar <- paste(idf1.dcast2$feature, idf1.dcast2$data)
#   colnames(idf1.dcast2)[3:4] <- paste(colnames(idf1.dcast2)[3:4],".FI", sep = "")
#   idf1.dcast <- merge(idf1.dcast, idf1.dcast2[,c("mvar","class")], by="mvar", all.x=TRUE)
#   idf1.dcast <- idf1.dcast[!is.na(idf1.dcast$class),]
#   valid.classes <- unique(idf1.dcast$class)
#   valid.comparisons <- combn(valid.classes, 2, simplify = FALSE)
#   print(valid.comparisons)
#   boxplots.featImp[[i]][["bininter"]][["mda"]] <- ggplot(idf1.dcast, aes(x=class, y=bininter.FI)) + 
#     geom_boxplot(outlier.shape = NA) + 
#     geom_point() + 
#     facet_grid(.~data, scales = "free_x", space = "free_x") + 
#     stat_compare_means(aes(group=class), 
#                        method = "wilcox.test", 
#                        comparisons = valid.comparisons,
#                        label = "p.signif", vjust = 1.5) + 
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none")
#   boxplots.featImp[[i]][["bininter"]][["indval"]] <- ggplot(idf1.dcast, aes(x=class, y=indval.FI)) + 
#     geom_boxplot(outlier.shape = NA) + 
#     geom_point() + 
#     facet_grid(.~data, scales = "free_x", space = "free_x") + 
#     stat_compare_means(aes(group=class), 
#                        method = "wilcox.test", 
#                        comparisons = valid.comparisons,
#                        label = "p.signif", vjust = 1.5) + 
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none")
#   ##indval-terinter
#   idf1 <- idf[idf$source %in% c("indval","terinter"),]
#   idf1.dcast <- dcast(data = idf1, formula = feature+data~source, value.var="featureImportance")
#   idf1.dcast2 <- dcast(data = idf1, formula = feature+data~source, value.var="IsIndSp")
#   idf1.dcast2$class <- ifelse((idf1.dcast2$terinter==1 & idf1.dcast2$indval==1),"common",
#                               ifelse((idf1.dcast2$terinter==1 & idf1.dcast2$indval==0),"terinterOnly",
#                                      ifelse(idf1.dcast2$terinter==0 & idf1.dcast2$indval==1,"indvalOnly",NA)))
#   idf1.dcast2 <- idf1.dcast2[!is.na(idf1.dcast2$class),]
#   ilist.terinter <- list()
#   for(j in unique(idf1.dcast2$class))
#   {
#     print(j)
#     for(m in unique(idf1.dcast2$data))
#     {
#       print(m)
#       ##Get the features 
#       jfeats <- idf1.dcast2$feature[idf1.dcast2$class %in% j & idf1.dcast2$data %in% m]
#       ##Correlate feature importance metrics if >1
#       if(length(jfeats)>1)
#       {
#         jfeats.corr <- cor.test(idf1.dcast$terinter[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m], 
#                                 idf1.dcast$indval[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m], 
#                                 method="spearman")
#         jfeats.corr.df <- data.frame(comparison="indval-terinter",
#                                      feature.class=j,
#                                      dataSource=m,
#                                      rho=jfeats.corr$estimate,
#                                      pval=jfeats.corr$p.value,
#                                      features=length(jfeats))
#         ilist.terinter[[paste(j,m)]] <- jfeats.corr.df
#       }
#     }
#   }
#   ilist.terinter <- do.call("rbind", ilist.terinter)
#   corr.list[[i]] <- rbind(ilist.terinter,ilist.bininter)
#   ##boxplots featureImportance by feature class
#   idf1.dcast$mvar <- paste(idf1.dcast$feature, idf1.dcast$data)
#   colnames(idf1.dcast)[3:4] <- paste(colnames(idf1.dcast)[3:4],".FI", sep = "")
#   idf1.dcast2$mvar <- paste(idf1.dcast2$feature, idf1.dcast2$data)
#   colnames(idf1.dcast2)[3:4] <- paste(colnames(idf1.dcast2)[3:4],".FI", sep = "")
#   idf1.dcast <- merge(idf1.dcast, idf1.dcast2[,c("mvar","class")], by="mvar", all.x=TRUE)
#   idf1.dcast <- idf1.dcast[!is.na(idf1.dcast$class),]
#   valid.classes <- unique(idf1.dcast$class)
#   valid.comparisons <- combn(valid.classes, 2, simplify = FALSE)
#   print(valid.comparisons)  # Debug
#   boxplots.featImp[[i]][["terinter"]][["mda"]] <- ggplot(idf1.dcast, aes(x=class, y=terinter.FI)) + 
#     geom_boxplot(outlier.shape = NA) + 
#     geom_point() + 
#     facet_grid(.~data, scales = "free_x", space = "free_x") + 
#     stat_compare_means(aes(group=class), 
#                        method = "wilcox", 
#                        comparisons = valid.comparisons,
#                        label = "p.signif", vjust = 1.5) + 
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none")
#   boxplots.featImp[[i]][["terinter"]][["indval"]] <- ggplot(idf1.dcast, aes(x=class, y=indval.FI)) + 
#     geom_boxplot(outlier.shape = NA) + 
#     geom_point() + 
#     facet_grid(.~data, scales = "free_x", space = "free_x") + 
#     stat_compare_means(aes(group=class), 
#                        method = "wilcox.test", 
#                        comparisons = valid.comparisons,
#                        label = "p.signif", vjust = 1.5) + 
#     theme_bw() + 
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           legend.position = "none")
# }

# Main update on boxplots to handle missing dataon bininterOnly class 
########################################################################################################################

boxplots.featImp <- list()
for(i in unique(mergedf.all$comparison)) {
  print(i)
  idf <- mergedf.all[mergedf.all$comparison %in% i,]
  
  ##indval-bininter
  idf1 <- idf[idf$source %in% c("indval","bininter"),]
  idf1.dcast <- dcast(data = idf1, formula = feature+data~source, value.var="featureImportance")
  idf1.dcast2 <- dcast(data = idf1, formula = feature+data~source, value.var="IsIndSp")
  idf1.dcast2$class <- ifelse((idf1.dcast2$bininter==1 & idf1.dcast2$indval==1),"common",
                              ifelse((idf1.dcast2$bininter==1 & idf1.dcast2$indval==0),"bininterOnly",
                                     ifelse(idf1.dcast2$bininter==0 & idf1.dcast2$indval==1,"indvalOnly",NA)))
  idf1.dcast2 <- idf1.dcast2[!is.na(idf1.dcast2$class),]
  
  ##boxplots featureImportance by feature class
  idf1.dcast$mvar <- paste(idf1.dcast$feature, idf1.dcast$data)
  colnames(idf1.dcast)[3:4] <- paste(colnames(idf1.dcast)[3:4],".FI", sep = "")
  idf1.dcast2$mvar <- paste(idf1.dcast2$feature, idf1.dcast2$data)
  colnames(idf1.dcast2)[3:4] <- paste(colnames(idf1.dcast2)[3:4],".FI", sep = "")
  idf1.dcast <- merge(idf1.dcast, idf1.dcast2[,c("mvar","class")], by="mvar", all.x=TRUE)
  idf1.dcast <- idf1.dcast[!is.na(idf1.dcast$class),]
  
  ilist.bininter <- list()
  for(j in unique(idf1.dcast2$class))
  {
    print(j)
    for(m in unique(idf1.dcast2$data))
    {
      print(m)
      ##Get the features 
      jfeats <- idf1.dcast2$feature[idf1.dcast2$class %in% j & idf1.dcast2$data %in% m]
      ##Correlate feature importance metrics if >1
      if(length(jfeats)>1)
      {
        jfeats.corr <- cor.test(idf1.dcast$bininter[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m], 
                                idf1.dcast$indval[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m], 
                                method="spearman")
        jfeats.corr.df <- data.frame(comparison="indval-bininter",
                                     feature.class=j,
                                     dataSource=m,
                                     rho=jfeats.corr$estimate,
                                     pval=jfeats.corr$p.value,
                                     features=length(jfeats))
        ilist.bininter[[paste(j,m)]] <- jfeats.corr.df
      }
    }
  }
  ilist.bininter <- do.call("rbind", ilist.bininter)
  # Create a function to generate comparisons for each facet
  generate_facet_comparisons <- function(facet_data) {
    valid_classes <- unique(facet_data$class)
    if(length(valid_classes) >= 2) {
      return(combn(valid_classes, 2, simplify = FALSE))
    } else {
      return(NULL)
    }
  }
  
  # Split data by facet and generate comparisons for each
  facet_data_list <- split(idf1.dcast, idf1.dcast$data)
  facet_comparisons <- lapply(facet_data_list, generate_facet_comparisons)
  
  # Create the plot with facet-specific comparisons
  boxplots.featImp[[i]][["bininter"]][["mda"]] <- ggplot(idf1.dcast, aes(x=class, y=bininter.FI)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point() + 
    facet_grid(.~data, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Add comparisons for each facet that has valid comparisons
  for(facet_name in names(facet_comparisons)) {
    if(!is.null(facet_comparisons[[facet_name]])) {
      boxplots.featImp[[i]][["bininter"]][["mda"]] <- boxplots.featImp[[i]][["bininter"]][["mda"]] +
        stat_compare_means(
          data = idf1.dcast[idf1.dcast$data == facet_name,],
          aes(group = class),
          method = "wilcox.test",
          comparisons = facet_comparisons[[facet_name]],
          label = "p.signif",
          vjust = 1.5
        )
    }
  }
  
  # Repeat for indval.FI plot
  boxplots.featImp[[i]][["bininter"]][["indval"]] <- ggplot(idf1.dcast, aes(x=class, y=indval.FI)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point() + 
    facet_grid(.~data, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  for(facet_name in names(facet_comparisons)) {
    if(!is.null(facet_comparisons[[facet_name]])) {
      boxplots.featImp[[i]][["bininter"]][["indval"]] <- boxplots.featImp[[i]][["bininter"]][["indval"]] +
        stat_compare_means(
          data = idf1.dcast[idf1.dcast$data == facet_name,],
          aes(group = class),
          method = "wilcox.test",
          comparisons = facet_comparisons[[facet_name]],
          label = "p.signif",
          vjust = 1.5
        )
    }
  }
  
  ##indval-terinter (similar modifications as above)
  idf1 <- idf[idf$source %in% c("indval","terinter"),]
  idf1.dcast <- dcast(data = idf1, formula = feature+data~source, value.var="featureImportance")
  idf1.dcast2 <- dcast(data = idf1, formula = feature+data~source, value.var="IsIndSp")
  idf1.dcast2$class <- ifelse((idf1.dcast2$terinter==1 & idf1.dcast2$indval==1),"common",
                              ifelse((idf1.dcast2$terinter==1 & idf1.dcast2$indval==0),"terinterOnly",
                                     ifelse(idf1.dcast2$terinter==0 & idf1.dcast2$indval==1,"indvalOnly",NA)))
  idf1.dcast2 <- idf1.dcast2[!is.na(idf1.dcast2$class),]
  ilist.terinter <- list()
  for(j in unique(idf1.dcast2$class))
  {
    print(j)
    for(m in unique(idf1.dcast2$data))
    {
      print(m)
      ##Get the features
      jfeats <- idf1.dcast2$feature[idf1.dcast2$class %in% j & idf1.dcast2$data %in% m]
      ##Correlate feature importance metrics if >1
      if(length(jfeats)>1)
      {
        jfeats.corr <- cor.test(idf1.dcast$terinter[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m],
                                idf1.dcast$indval[idf1.dcast$feature %in% jfeats & idf1.dcast$data %in% m],
                                method="spearman")
        jfeats.corr.df <- data.frame(comparison="indval-terinter",
                                     feature.class=j,
                                     dataSource=m,
                                     rho=jfeats.corr$estimate,
                                     pval=jfeats.corr$p.value,
                                     features=length(jfeats))
        ilist.terinter[[paste(j,m)]] <- jfeats.corr.df
      }
    }
  }
  ilist.terinter <- do.call("rbind", ilist.terinter)
  corr.list[[i]] <- rbind(ilist.terinter,ilist.bininter)
  ##boxplots featureImportance by feature class
  idf1.dcast$mvar <- paste(idf1.dcast$feature, idf1.dcast$data)
  colnames(idf1.dcast)[3:4] <- paste(colnames(idf1.dcast)[3:4],".FI", sep = "")
  idf1.dcast2$mvar <- paste(idf1.dcast2$feature, idf1.dcast2$data)
  colnames(idf1.dcast2)[3:4] <- paste(colnames(idf1.dcast2)[3:4],".FI", sep = "")
  idf1.dcast <- merge(idf1.dcast, idf1.dcast2[,c("mvar","class")], by="mvar", all.x=TRUE)
  idf1.dcast <- idf1.dcast[!is.na(idf1.dcast$class),]
  
  # Split data by facet and generate comparisons for each
  facet_data_list <- split(idf1.dcast, idf1.dcast$data)
  facet_comparisons <- lapply(facet_data_list, generate_facet_comparisons)
  
  # Create terinter plots with facet-specific comparisons
  boxplots.featImp[[i]][["terinter"]][["mda"]] <- ggplot(idf1.dcast, aes(x=class, y=terinter.FI)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point() + 
    facet_grid(.~data, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  for(facet_name in names(facet_comparisons)) {
    if(!is.null(facet_comparisons[[facet_name]])) {
      boxplots.featImp[[i]][["terinter"]][["mda"]] <- boxplots.featImp[[i]][["terinter"]][["mda"]] +
        stat_compare_means(
          data = idf1.dcast[idf1.dcast$data == facet_name,],
          aes(group = class),
          method = "wilcox.test",
          comparisons = facet_comparisons[[facet_name]],
          label = "p.signif",
          vjust = 1.5
        )
    }
  }
  
  boxplots.featImp[[i]][["terinter"]][["indval"]] <- ggplot(idf1.dcast, aes(x=class, y=indval.FI)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point() + 
    facet_grid(.~data, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  for(facet_name in names(facet_comparisons)) {
    if(!is.null(facet_comparisons[[facet_name]])) {
      boxplots.featImp[[i]][["terinter"]][["indval"]] <- boxplots.featImp[[i]][["terinter"]][["indval"]] +
        stat_compare_means(
          data = idf1.dcast[idf1.dcast$data == facet_name,],
          aes(group = class),
          method = "wilcox.test",
          comparisons = facet_comparisons[[facet_name]],
          label = "p.signif",
          vjust = 1.5
        )
    }
  }
}

#######################################################################################################################

corr.list.inoff <- corr.list$INSHORE_OFFSHORE
corr.list.inoff$pval.cat <- ifelse(corr.list.inoff$pval<0.05,"*","")

corr.plots <- list()
for(i in names(corr.list))
{
  print(i)
  idf <- corr.list[[i]]
  idf$pval.cat <- ifelse(idf$pval<0.05,"*","")
  idf$comparison <- gsub("-","\nvs.\n", idf$comparison)
  corr.plots[[i]] <- ggplot(idf, aes(x=feature.class, y=dataSource)) + 
    geom_point(aes(size=features, colour = rho)) + 
    geom_text(aes(label=pval.cat), colour="white") + 
    scale_color_gradient2() + 
    ggtitle(i) + 
    xlab("Indicator Species set") + 
    ylab("Source data") + 
    labs(colour="Rho feature Importance") +
    facet_grid(.~comparison, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
do.call("grid.arrange", c(corr.plots, ncol=2, nrow=2))

grid.arrange(boxplots.featImp$INSHORE_OFFSHORE$bininter$mda,
             boxplots.featImp$INSHORE_OFFSHORE$bininter$indval,
             boxplots.featImp$INSHORE_OFFSHORE$terinter$mda,
             boxplots.featImp$INSHORE_OFFSHORE$terinter$indval, ncol=2, nrow=2)

# hlay <- rbind(c(1,1,2,2,3,3),
#               c(4,4,5,5,6,7))
# pdf(file="/data/projects/aime/analyses/bruvs/Figure3/Figure3.pdf", h=8, w=12)
# grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
#              plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$maxN + ggtitle("B"),
#              plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$`pres/abs` + ggtitle("C"),
#              prevplots$INSHORE_OFFSHORE + ggtitle("D") + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#              corr.plots$INSHORE_OFFSHORE + ggtitle("E") + labs(colour="Rho feature Importance"),
#              boxplots.featImp$INSHORE_OFFSHORE$bininter$mda + ggtitle("F"),
#              boxplots.featImp$INSHORE_OFFSHORE$terinter$mda + ggtitle("G"), layout_matrix=hlay)
# dev.off()
# 
# hlay <- rbind(c(1,1,2,3,4,4),
#               c(5,5,6,6,7,8))
# pdf(file="/data/projects/aime/analyses/bruvs/Figure3/Figure3.pdf", h=8, w=12)
# grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
#              plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$maxN + ggtitle("B"),
#              plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$`pres/abs` + ggtitle("C"),
#              adonis.plots$INSHORE_OFFSHORE + ggtitle("D"),
#              prevplots$INSHORE_OFFSHORE + ggtitle("E") + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#              corr.plots$INSHORE_OFFSHORE + ggtitle("F") + labs(colour="Rho feature Importance"),
#              boxplots.featImp$INSHORE_OFFSHORE$bininter$mda + ggtitle("G"),
#              boxplots.featImp$INSHORE_OFFSHORE$terinter$mda + ggtitle("H"), layout_matrix=hlay)
# dev.off()
# 
# 
# hlay <- rbind(c(1,1,2,3,4,4),
#               c(1,1,NA,NA,4,4),
#               c(5,5,6,6,7,8),
#               c(5,5,6,6,7,8))
# pdf(file="/data/projects/aime/analyses/bruvs/Figure3/Figure3_A.pdf", h=8, w=12)
# grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
#              plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$maxN + ggtitle("B"),
#              plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$`pres/abs` + ggtitle("C"),
#              adonis.plots$INSHORE_OFFSHORE + ggtitle("D"),
#              prevplots$INSHORE_OFFSHORE + ggtitle("E") + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#              corr.plots$INSHORE_OFFSHORE + ggtitle("F") + labs(colour="Rho feature Importance"),
#              boxplots.featImp$INSHORE_OFFSHORE$bininter$mda + ggtitle("G"),
#              boxplots.featImp$INSHORE_OFFSHORE$terinter$mda + ggtitle("H"), layout_matrix=hlay)
# dev.off()


hlay <- rbind(c(1,1,1,2,3,4,8,8,8),
              c(1,1,1,5,6,7,8,8,8),
              c(9,9,10,10,10,11,11,12,12),
              c(9,9,10,10,10,11,11,12,12))
pdf(file="/data/projects/aime/analyses/bruvs/Figure3/Figure3_Best.pdf", h=10, w=15)
grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$maxN + ggtitle("B"),
             nullGrob(),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$`pres/abs` + ggtitle("C"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byData$bininter + ggtitle("D"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byData$indval + ggtitle("E"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byData$terinter + ggtitle("F"),
             adonis.plots$INSHORE_OFFSHORE + ggtitle("G"),
             prevplots$INSHORE_OFFSHORE + ggtitle("H") ,
             corr.plots$INSHORE_OFFSHORE + ggtitle("I"),
             boxplots.featImp$INSHORE_OFFSHORE$bininter$mda + ggtitle("J"),
             boxplots.featImp$INSHORE_OFFSHORE$terinter$mda + ggtitle("K"), layout_matrix=hlay)
dev.off()

# library(patchwork)
# hlay <- 
#   "AABCDHH
#  AAEFGHH
#  IIJJJKL"
# pdf(file="/data/projects/aime/analyses/bruvs/Figure3/Figure3_Best.pdf", h=8, w=12)
# plots.list$INSHORE_OFFSHORE$barplot +
#   plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$maxN + plot_spacer() +
#   plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$`pres/abs` +
#   plots.list$INSHORE_OFFSHORE$vennDiagrams_byData$bininter +
#   plots.list$INSHORE_OFFSHORE$vennDiagrams_byData$indval +
#   plots.list$INSHORE_OFFSHORE$vennDiagrams_byData$terinter + 
#   adonis.plots$INSHORE_OFFSHORE + 
#   prevplots$INSHORE_OFFSHORE + 
#   corr.plots$INSHORE_OFFSHORE +
#   boxplots.featImp$INSHORE_OFFSHORE$bininter$mda + 
#   boxplots.featImp$INSHORE_OFFSHORE$terinter$mda + plot_layout(design = hlay) + plot_annotation("A")
# dev.off()

##Get the stats for the pres/abs data (just 2 levels in class, ggpubr with facets don't plot it)
df <- boxplots.featImp$INSHORE_OFFSHORE$bininter$mda$data
df <- df[df$data %in% "pres/abs",]
wilcox.test(df$bininter.FI~df$class)
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  df$bininter.FI by df$class
# W = 109, p-value = 0.002745
# alternative hypothesis: true location shift is not equal to 0

##From ggpubr::stat_compare_means, significativity level=**=p<0.01
#Added to the pdf plot with inkscape

#------------------

plots.list.barplot.all <- lapply(plots.list, function(x){x[["barplot"]][["data"]]})
plots.list.barplot.all <- do.call("rbind", plots.list.barplot.all)
plots.list.barplot.all.plot <- ggplot(plots.list.barplot.all[!plots.list.barplot.all$Var3 %in% "INSHORE_OFFSHORE",], aes(x=Var2, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_text(aes(label=Freq), vjust=0, position=position_dodge(.9)) + 
  scale_fill_manual(values = colours.method) + 
  ylab("Indicator Species") + 
  xlab("Data Source") + 
  facet_grid(.~Var3) + 
  labs(fill="Method") + 
  theme_bw()


prevplots.all <- lapply(prevplots, function(x){x[["data"]]})
prevplots.all <- do.call("rbind", prevplots.all)
prevplots.all.plot <- ggplot(prevplots.all[!prevplots.all$comparison %in% "INSHORE_OFFSHORE",], aes(x=source, y=prevalence)) + 
  geom_boxplot(aes(fill=IsIndSp), outlier.shape = NA, alpha=.5) + 
  geom_point(aes(colour=IsIndSp), position = position_dodge(width = .8)) + 
  facet_grid(comparison~data) + 
  scale_color_manual(values = c("black", "grey")) +
  scale_fill_manual(values = c("black", "grey")) +
  stat_compare_means(aes(group=IsIndSp), method = "wilcox", label = "p.signif", vjust =1.5) + 
  ylab("Species prevalence") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



corr.plots.all <- lapply(corr.plots, function(x){x[["data"]]})
for(i in names(corr.plots.all))
{
  corr.plots.all[[i]][,"SiteComparison"] <- i
}

corr.plots.all <- do.call("rbind",corr.plots.all)
corr.plots.all.plot <- ggplot(corr.plots.all[!corr.plots.all$SiteComparison %in% "INSHORE_OFFSHORE",], aes(x=feature.class, y=dataSource)) + 
  geom_point(aes(size=features, colour = rho)) + 
  geom_text(aes(label=pval.cat), colour="white") + 
  scale_color_gradient2() + 
  ggtitle(i) + 
  xlab("Indicator Species set") + 
  ylab("Source data") + 
  facet_grid(SiteComparison~comparison, scales = "free_x", space = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 360))



boxplots.featImp.allBin.mda <- lapply(boxplots.featImp, function(x){x[["bininter"]][["mda"]][["data"]]})
for(i in names(boxplots.featImp.allBin.mda))
{
  boxplots.featImp.allBin.mda[[i]][,"comparison"] <- i
}
boxplots.featImp.allBin.mda <- do.call("rbind", boxplots.featImp.allBin.mda)
boxplots.featImp.allBin.mda.plot <- ggplot(boxplots.featImp.allBin.mda[!boxplots.featImp.allBin.mda$comparison %in% "INSHORE_OFFSHORE",], aes(x=class, y=bininter.FI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  facet_grid(comparison~data, scales = "free_x", space = "free_x") + 
  stat_compare_means(aes(group=class), 
                     method = "wilcox", 
                     comparisons = combn(x = unique(boxplots.featImp.allBin.mda$class), m = 2, simplify = FALSE),
                     label = "p.signif", 
                     vjust = 1.5
                     ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


boxplots.featImp.allTer.mda <- lapply(boxplots.featImp, function(x){x[["terinter"]][["mda"]][["data"]]})
for(i in names(boxplots.featImp.allTer.mda))
{
  boxplots.featImp.allTer.mda[[i]][,"comparison"] <- i
}
boxplots.featImp.allTer.mda <- do.call("rbind", boxplots.featImp.allTer.mda)
boxplots.featImp.allTer.mda.plot <- ggplot(boxplots.featImp.allTer.mda[!boxplots.featImp.allTer.mda$comparison %in% "INSHORE_OFFSHORE",], aes(x=class, y=terinter.FI)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() + 
  facet_grid(comparison~data, scales = "free_x", space = "free_x") + 
  stat_compare_means(aes(group=class), 
                     method = "wilcox", 
                     comparisons = combn(x = unique(boxplots.featImp.allTer.mda$class), m = 2, simplify = FALSE),
                     label = "p.signif", 
                     vjust = 1.5
                     ) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


##integral visualization all adonis results
adonisplot <- ggplot(adonis.alldf[!adonis.alldf$comparison %in% "INSHORE_OFFSHORE",], aes(x=source, y=R2)) +
  geom_bar(aes(fill=features), stat="identity") +
  facet_grid(data~comparison) +
  geom_text(aes(label=pval.cat)) + 
  scale_fill_viridis() + 
  coord_flip() + 
  theme_bw()



# pdf(file="/data/projects/aime/analyses/bruvs/Figure3/FigureS2.pdf", h=15, w=10)
# hlay <- rbind(c(1,1,1),
#               c(2,2,NA),
#               c(3,4,5),
#               c(3,4,5))
# grid.arrange(plots.list.barplot.all.plot + ggtitle("A"),
#              corr.plots.all.plot + ggtitle("B"),
#              prevplots.all.plot + ggtitle("C"),
#              boxplots.featImp.allBin.mda.plot + ylab("Mean Decrease Accuracy; Bininter") + ggtitle("D"),
#              boxplots.featImp.allTer.mda.plot + ylab("Mean Decrease Accuracy; Terinter") + ggtitle("E"),
#              layout_matrix=hlay)
# dev.off()


pdf(file="/data/projects/aime/analyses/bruvs/Figure3/FigureS2_Best.pdf", h=15, w=10)
hlay <- rbind(c(1,1,1,1,1,1),
              c(2,2,2,3,3,3),
              c(4,4,5,5,6,6),
              c(4,4,5,5,6,6))
grid.arrange(plots.list.barplot.all.plot + ggtitle("A"),
             adonisplot + ggtitle("B") + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom"),
             corr.plots.all.plot + ggtitle("C"),
             prevplots.all.plot + ggtitle("D"),
             boxplots.featImp.allBin.mda.plot + ylab("Mean Decrease Accuracy; Bininter") + ggtitle("E"),
             boxplots.featImp.allTer.mda.plot + ylab("Mean Decrease Accuracy; Terinter") + ggtitle("F"),
             layout_matrix=hlay)
dev.off()

