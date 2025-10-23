library(ggplot2) # for plots
library(gridExtra) # for multiple plots
library(patchwork) # to display multiples graphs in the same figure
library(readr)
library(ggpubr)

dir.create("/data/projects/aime/analyses/bruvs/FigureS1/")

objlist <- list()

# load indval results
x <- load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/Indval_output_data/Indval_all_analyses_overall_data_prev_0.Rda")
for(i in x)
{
  objlist[["indval"]][[i]] <- get(i)
}
rm(list = x)

# load bininter results
x <- load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/bininter_output_data/bininter_Predomics_all_analyses_overall_data_prev_0.Rda")
for(i in x)
{
  objlist[["bininter"]][[i]] <- get(i)
}
rm(list = x)

# load terinter results
x <- load(file="/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/terinter_output_data/terinter_Predomics_all_analyses_overall_data_prev_0.Rda")
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
  plots.list[[i]][["vennDiagrams"]] <- lapply(i.listVenn, function(x){ggvenn(x,fill_color = c("#648C67","#8195B2","#CC3837"), stroke_size = 0.5, set_name_size = 4,show_percentage = FALSE, stroke_alpha = 1) + 
      theme(
        plot.title = element_text(hjust = 0.5),
        # legend.title = element_text(color = "black", size = 9, margin = margin(0, 10, 0, 0), hjust = 0.5),
        legend.position = "none"
      )  + 
      scale_x_continuous(expand = expansion(mult = 0.02))})
  
  i.listVenn2 <- list()
  for(j in unique(idf$source))
  {
    i.listVenn2[[j]] <- list(maxN=idf$feature[idf$source %in% j & idf$data %in% "maxN"],
                            presabs=idf$feature[idf$source %in% j & idf$data %in% "pres/abs"])
  }
  plots.list[[i]][["vennDiagrams_byMethod"]] <- lapply(i.listVenn2, function(x){ggvenn(x, stroke_size = 0.5, set_name_size = 4,show_percentage = FALSE, stroke_alpha = 1) + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        # legend.title = element_text(color = "black", size = 9, margin = margin(0, 10, 0, 0), hjust = 0.5),
        legend.position = "none"
      )  + 
      scale_x_continuous(expand = expansion(mult = 0.2))})
}

##visu the results of venn comparing abundance vs. presence absence for each method
grid.arrange(plots.list$INSHORE_OFFSHORE$barplot,
             plots.list$INSHORE_OFFSHORE$vennDiagrams$maxN,
             plots.list$INSHORE_OFFSHORE$vennDiagrams$`pres/abs`, ncol=3)

do.call("grid.arrange",c(plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod, ncol=3))

x <- list.files(path = "/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/bininter_output_data/", pattern = "*_prev_*", full.names = TRUE)
bininterObj <- list()
for(i in x)
{
  print(i)
  ilab <- gsub(".*_prev|\\.Rda","", basename(i))
  idata <- load(i)
  for(o in idata)
  {
    bininterObj[[ilab]][[o]] <- get(o)
  }
  rm(list=idata)
}

x <- list.files(path = "/data/projects/aime/analyses/metabarcoding/4.data_video_analysis/1.Article_figures/fig2/fig2-B-C-E-F/terinter_output_data/", pattern = "*_prev_*", full.names = TRUE)
x <- x[-1]

terinterObj <- list()
for(i in x)
{
  print(i)
  ilab <- gsub(".*_prev|\\.Rda","", basename(i))
  idata <- load(i)
  for(o in idata)
  {
    terinterObj[[ilab]][[o]] <- get(o)
  }
  rm(list=idata)
}

##bininter objects processing
tframes <- c("empirical.auc","generalization.auc","empirical.acc","generalization.acc")
bininter.metricsdf <- list()
for(i in names(bininterObj))
{
  print(i)
  for(j in c("predout.bin","predout.maxn"))
  {
    print(names(bininterObj[[i]][[j]]))
    for(c in names(bininterObj[[i]][[j]]))
    {
      print(names(bininterObj[[i]][[j]][[c]]))
      ijc.dfs <- bininterObj[[i]][[j]][[c]][["fit"]][["crossVal"]][["scores"]][tframes]
      ijc.dfs <- lapply(ijc.dfs, function(x){x[,"k"] <- rownames(x) ; return(melt(x))})
      for(l in names(ijc.dfs))
      {
        ijc.dfs[[l]][,"metric"] <- l
      }
      ijc.dfs <- do.call("rbind", ijc.dfs)
      ijc.dfs$prevalence <- i
      ijc.dfs$datasource <- j
      ijc.dfs$comparison <- c
      ijc.dfs$algo <- "bininter"
      bininter.metricsdf[[paste(i,j,c)]] <- ijc.dfs
    }
  }
}
bininter.metricsdf <- do.call("rbind",bininter.metricsdf)

##bininter objects processing
tframes <- c("empirical.auc","generalization.auc","empirical.acc","generalization.acc")
terinter.metricsdf <- list()
for(i in names(terinterObj))
{
  print(i)
  for(j in c("predout.bin","predout.maxn"))
  {
    print(names(terinterObj[[i]][[j]]))
    for(c in names(terinterObj[[i]][[j]]))
    {
      print(names(terinterObj[[i]][[j]][[c]]))
      ijc.dfs <- terinterObj[[i]][[j]][[c]][["fit"]][["crossVal"]][["scores"]][tframes]
      ijc.dfs <- lapply(ijc.dfs, function(x){x[,"k"] <- rownames(x) ; return(melt(x))})
      for(l in names(ijc.dfs))
      {
        ijc.dfs[[l]][,"metric"] <- l
      }
      ijc.dfs <- do.call("rbind", ijc.dfs)
      ijc.dfs$prevalence <- i
      ijc.dfs$datasource <- j
      ijc.dfs$comparison <- c
      ijc.dfs$algo <- "terinter"
      terinter.metricsdf[[paste(i,j,c)]] <- ijc.dfs
    }
  }
}
terinter.metricsdf <- do.call("rbind",terinter.metricsdf)

metricsdfs.pred <- rbind(terinter.metricsdf,bininter.metricsdf)
metricsdfs.pred$datasource <- factor(metricsdfs.pred$datasource, levels=c("predout.maxn","predout.bin"), labels=c("maxN","presAbs"))
metricsdfs.pred$prevalence <- factor(metricsdfs.pred$prevalence, levels=c("_0","_10"), labels = c("all\nspecies","species >10%\nprevalence"))

vplots <- list()
for(i in unique(metricsdfs.pred$comparison))
{
  print(i)
  idf <- metricsdfs.pred[metricsdfs.pred$comparison %in% i & metricsdfs.pred$metric %in% c("empirical.auc","generalization.auc"),]
  idf$metric <- factor(idf$metric, levels=c("empirical.auc","generalization.auc"), labels=c("Empirical","Generalization"))
  vplots[[i]] <- ggplot(idf, aes(x=metric, y=value)) +
    geom_violin(aes(fill=prevalence), trim = TRUE, alpha = 0.5, position = position_dodge(width = 0.9), scale = "width") +
    geom_boxplot(aes(fill=prevalence), width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9), alpha = 0.7) +
    ylim(0,1.1) +
    facet_grid(algo~datasource) + 
    ylab("AUC") + 
    stat_compare_means(aes(group=prevalence), method = "wilcox", label = "p.signif", vjust=1) +
    theme_bw()
}

hlay <- rbind(c(1,2,3),
              c(1,NA,NA),
              c(4,4,4),
              c(4,4,4))
pdf(file="/data/projects/aime/analyses/bruvs/FigureS1/FigS1.pdf", h=10, w=10)
grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams$maxN + ggtitle("B"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams$`pres/abs` + ggtitle("C"), 
             vplots$INSHORE_OFFSHORE + ggtitle("D"), layout_matrix=hlay)
dev.off()

hlay <- rbind(c(1,2,3,4),
              c(1,NA,NA,NA),
              c(5,5,5,5),
              c(5,5,5,5))
pdf(file="/data/projects/aime/analyses/bruvs/FigureS1/FigS1_v2.pdf", h=10, w=10)
grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$indval + ggtitle("B"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$bininter + ggtitle("C"), 
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$terinter + ggtitle("D"), 
             vplots$INSHORE_OFFSHORE + ggtitle("E"), layout_matrix=hlay)
dev.off()

hlay <- rbind(c(1,1,2,3,4),
              c(1,1,5,6,7),
              c(8,8,8,8,8),
              c(8,8,8,8,8))

pdf(file="/data/projects/aime/analyses/bruvs/FigureS1/FigureS1_Best.pdf", h=10, w=10)
grid.arrange(plots.list$INSHORE_OFFSHORE$barplot + ggtitle("A"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams$maxN + ggtitle("B"),
             nullGrob(),
             plots.list$INSHORE_OFFSHORE$vennDiagrams$`pres/abs` + ggtitle("C"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$bininter + ggtitle("D"), 
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$indval + ggtitle("E"),
             plots.list$INSHORE_OFFSHORE$vennDiagrams_byMethod$terinter + ggtitle("F"),
             vplots$INSHORE_OFFSHORE + ggtitle("G"), layout_matrix=hlay)
dev.off()
