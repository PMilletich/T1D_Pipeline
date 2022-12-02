library(effectsize)
library(ggplot2)
library(ggpubr)

###
#Picrust output
###########
Sample_data = read.csv("./Matched_Samples.csv")
Sample_data_1 = Sample_data[,c("ID", "Autoimmune_2_groups")]
Pathway = read.table("Picrust_totalPop_TotalAbun/pathways_out/path_abun_unstrat_descrip.tsv")
colnames(Pathway) = Pathway[1,]
Pathway = Pathway[-1,]
rownames(Pathway) = Pathway$pathway
Pathway_Descript = Pathway[,c("pathway", "description")]
colnames(Pathway_Descript) = c("Pathway", "description")
Pathway$`pathway` = NULL
Pathway$description = NULL
Pathway_1 = Pathway
for (i in colnames(Pathway_1)) {
  new_col = Pathway_1[,i]
  new_col_1 = as.numeric(new_col)/(sum(as.numeric(new_col)))
  Pathway_1[,i] = new_col_1
}
Pathway_Samples = data.frame(t(Pathway_1))

Pathway_Names = colnames(Pathway_Samples)
Pathway_Samples$ID = gsub("X", "", rownames(Pathway_Samples))
Pathway_Samples = merge(Sample_data_1, Pathway_Samples)

current_pathway = "PWY.3781"
graph_data = data.frame()
cohens_data = data.frame()
fdr_data = data.frame()
fdr_data_t = data.frame()
fdr_data_w = data.frame()
for (current_pathway in Pathway_Names) {
  current_data = Pathway_Samples[,c("Autoimmune_2_groups", current_pathway)]
  colnames(current_data)  = c("Autoimmune_2_groups", "Pathway")
  current_data$Pathway = as.numeric(current_data$Pathway)
  
  Control = subset(current_data, current_data$Autoimmune_2_groups == "Control")$Pathway
  T1D = subset(current_data, current_data$Autoimmune_2_groups != "Control")$Pathway

  Normalcy = shapiro.test(current_data$Pathway)
  Normalcy = Normalcy$p.value
  
  Wilcox_output = wilcox.test(Pathway ~ Autoimmune_2_groups, data = current_data)
  Ttest_output = t.test(T1D, Control)
  T_Pvalue = Ttest_output$p.value
  W_Pvalue = Wilcox_output$p.value
  
  Effect = round(cohens_d(Pathway ~ Autoimmune_2_groups, data = current_data)$Cohens_d,3)

  if (Normalcy > 0.05) {
    if (T_Pvalue <= 0.05) {
      print(paste("Ttest;", current_pathway, Effect))
      new_data = data.frame("Abundance"= Control, "Group" = "Control", "Pathway" = current_pathway)
      graph_data = rbind(graph_data, new_data)
      new_data = data.frame("Abundance"= T1D, "Group" = "Type 1 Diabetes", "Pathway" = current_pathway)
      graph_data = rbind(graph_data, new_data)
      
      new_data = data.frame("EffectSize"= Effect, "Pathway" = current_pathway)
      cohens_data = rbind(cohens_data, new_data)
      
    }
    current_data = data.frame("Pathway" = current_pathway, "Pvalue" = T_Pvalue)
    fdr_data_t = rbind(fdr_data_t, current_data)
  } else {
    if (W_Pvalue <= 0.05) {
      print(paste("Wilcox;", current_pathway, Effect))
      new_data = data.frame("Abundance"= Control, "Group" = "Control", "Pathway" = current_pathway)
      graph_data = rbind(graph_data, new_data)
      new_data = data.frame("Abundance"= T1D, "Group" = "Type 1 Diabetes", "Pathway" = current_pathway)
      graph_data = rbind(graph_data, new_data)
      
      new_data = data.frame("EffectSize"= Effect, "Pathway" = current_pathway)
      cohens_data = rbind(cohens_data, new_data)
    }
    current_data = data.frame("Pathway" = current_pathway, "Pvalue" = W_Pvalue)
    fdr_data_w = rbind(fdr_data_w, current_data)
  }
}

fdr_data$fdr = p.adjust(fdr_data$Pvalue, "fdr")
fdr_data_t$fdr = p.adjust(fdr_data_t$Pvalue, "fdr")
fdr_data_w$fdr = p.adjust(fdr_data_w$Pvalue, "fdr")

graph_data$Pathway = gsub("\\.", "-", graph_data$Pathway)
graph_data_1 = merge(Pathway_Descript, graph_data)

graph_data_1$description= gsub("acetyl-CoA_fermentation_to_butanoate_II", 
     "acetyl-CoA fermentation to \nbutanoate II", graph_data_1$description)                       
graph_data_1$description= gsub("pyruvate_fermentation_to_acetone", 
                               "pyruvate fermentation to acetone", graph_data_1$description) 
graph_data_1$description= gsub("yrinate_a,c-diamide_biosynthesis_I_", 
                               "yrinate a,c-diamide biosynthesis I\n",
                               graph_data_1$description) 
graph_data_1$description= gsub("nitrate_reduction_VI_",
                               "nitrate reduction VI ", graph_data_1$description) 
graph_data_1$description = gsub("_", " ", graph_data_1$description)


Abundance = ggplot(graph_data_1, aes(x = description, y = Abundance, fill = Group)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 1) + ylab("Relative Abundance") + xlab("") +
  stat_compare_means( hide.ns = T, method = "wilcox", size = 3,label = "p.format",
                      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
                                         symbols = c("****", "***", "**", "*", ".", "ns")))+ 
  ylim(c(0,0.008)) + 
  theme(axis.text.x = element_text(size = 6),
        legend.position="top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(breaks = c("Control", "Type 1 Diabetes" ),
                    values = c("#888888", "#6699CC")) + 
  annotate("text", x = 1, y = 0.0078, label = "Cohen's d = 0.511", size = 3)+ 
  annotate("text", x = 2, y = 0.0078, label = "Cohen's d = 0.447", size = 3)+ 
  annotate("text", x = 3, y = 0.0078, label = "Cohen's d = 0.552", size = 3)+ 
  annotate("text", x = 4, y = 0.0078, label = "Cohen's d = 0.417", size = 3) + 
  annotate("text", x = 1, y = 0.0075, label = "padj: ns", size = 3)+ 
  annotate("text", x = 2, y = 0.0075, label = "padj: ns", size = 3)+ 
  annotate("text", x = 3, y = 0.0075, label = "padj: ns", size = 3)+ 
  annotate("text", x = 4, y = 0.0075, label = "padj: ns", size = 3); Abundance


cohens_data$Pathway = gsub("\\.", "-", cohens_data$Pathway)
cohens_data = merge(Pathway_Descript, cohens_data)
cohens_data$description= gsub("acetyl-CoA_fermentation_to_butanoate_II", 
                               "acetyl-CoA fermentation to \nbutanoate II", cohens_data$description)                       
cohens_data$description= gsub("pyruvate_fermentation_to_acetone", 
                               "pyruvate fermentation to acetone", cohens_data$description) 
cohens_data$description= gsub("yrinate_a,c-diamide_biosynthesis_I_", 
                               "yrinate a,c-diamide biosynthesis I\n",
                              cohens_data$description) 
cohens_data$description= gsub("nitrate_reduction_VI_",
                               "nitrate reduction VI ", cohens_data$description) 

jpeg("Pathways.jpeg", res = 600, height = 4000, width = 4000)
Abundance
dev.off()
