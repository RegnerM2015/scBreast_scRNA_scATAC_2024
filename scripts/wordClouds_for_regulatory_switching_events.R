################################################################################
# Matt Regner
# Franco Lab
# Plot word clouds for regulatory switching events

# Run in RStudio not docker container

library(RColorBrewer)
library(wordcloud)

# Make output folder
if(dir.exists("wordCloud_Plots")){
  print("Directory wordCloud_Plots already exists!")
}else{
  dir.create("wordCloud_Plots")
}


# Make word cloud for Basal +/- upregulated genes
p2g <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/p2gs_interaction.rds")

# Read in genes upregulated in cancer v. normal
up.genes <- readRDS("./Basal_Cohort_Results-updates-Revised-updated/up_genes.rds")
p2g <- p2g[p2g$DiffClass == "+/-" & p2g$OLS_model_geneName %in% up.genes$gene,]


df <- data.frame(table(p2g$LME_cell_type_1_model_geneName))

pdf(paste0("wordCloud_Plots/word_cloud_Basal_-_normal_+_cancer-upgenes-",nrow(df),"unique_genes.pdf"))
wordcloud(words = df$Var1, freq = df$Freq, min.freq = 1,
          max.words=500, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))
dev.off()


# Make word cloud for Luminal +/- upregulated genes
p2g <- readRDS("./Luminal_Cohort_Results-updates-Revised-updated/p2gs_interaction.rds")

# Read in genes upregulated in cancer v. normal
up.genes <- readRDS("./Luminal_Cohort_Results-updates-Revised-updated/up_genes.rds")
p2g <- p2g[p2g$DiffClass == "+/-" & p2g$OLS_model_geneName %in% up.genes$gene,]


df <- data.frame(table(p2g$LME_cell_type_1_model_geneName))

pdf(paste0("wordCloud_Plots/word_cloud_Luminal_-_normal_+_cancer-upgenes-",nrow(df),"unique_genes.pdf"))
wordcloud(words = df$Var1, freq = df$Freq, min.freq = 1,
          max.words=500, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))
dev.off()
