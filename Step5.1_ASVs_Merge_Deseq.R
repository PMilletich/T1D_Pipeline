library(ggplot2)
library(ggpubr)
setwd("~/Desktop/T1D_Demultiplexed/")

taxa_identifiers = read.csv("Taxa_Identifiers_v2.csv")
taxa_identifiers_1 = data.frame()
current_glom = "Identifier"
current_method= "Reads/g"
iteration_Threshold = 50

for (current_method in c("Reads/g", "Relative Abundance")) {
  current_reads = read.csv(paste("./CSV_files/Deseq_Reads_", current_glom, "_2.csv", sep = ""))
  reads_Count = data.frame(table(current_reads$Taxa_Name))
  reads_Count = subset(reads_Count, reads_Count$Freq >= iteration_Threshold)
  Reads_Taxa_list = unique(reads_Count$Var1)
  Reads_Taxa_list = subset(Reads_Taxa_list, grepl("NA.", Reads_Taxa_list) == F)
  
  current_rel = read.csv(paste("./CSV_files/Deseq_RelAbun_", current_glom, "_2.csv", sep = ""))
  RelAbun_Count = data.frame(table(current_rel$Taxa_Name))
  RelAbun_Count = subset(RelAbun_Count, RelAbun_Count$Freq >= iteration_Threshold)
  RelAbun_Taxa_list = unique(RelAbun_Count$Var1)
  RelAbun_Taxa_list = subset(RelAbun_Taxa_list, grepl("NA.", RelAbun_Taxa_list) == F)
  
  Taxa_list = subset(Reads_Taxa_list, Reads_Taxa_list %in% RelAbun_Taxa_list)
  
  if (current_method == "Reads/g") {
    current_data = current_reads
  } else {
    current_data = current_rel
  }

  current_Count = data.frame(table(current_data$Taxa_Name))
  current_Count = subset(current_Count, current_Count$Freq >= iteration_Threshold)
  Taxa_list_OG = unique(current_Count$Var1)
  Taxa_list_OG = subset(Taxa_list_OG, grepl("NA_NA", Taxa_list_OG) == F)
  
  print(length(Taxa_list_OG))
  print(length(unique(current_data$Taxa_Name)))
  current_data = subset(current_data, current_data$Taxa_Name %in% Taxa_list_OG)
  print(length(unique(current_data$Taxa_Name)))
  current_data = subset(current_data, current_data$Taxa_Name %in% Taxa_list)
  print(length(unique(current_data$Taxa_Name)))
  
  df_Count = data.frame(table(current_data$Taxa_Name))
  colnames(df_Count) = c("Taxa_Name", "Iterations")
  current_data_Condensed = data.frame()
  
  current_data = merge(df_Count, current_data)
  
  current_data$Taxa = NULL
  
  current_Taxa = unique(current_data$Taxa_Name)[1]
  for (current_Taxa in unique(current_data$Taxa_Name)) {
    current_subset = subset(current_data, current_data$Taxa_Name == current_Taxa)
    
    current_LFC = mean(current_subset$log2FoldChange)
    current_LFC_SD = sd(current_subset$log2FoldChange)
    
    current_AI_prev = mean(current_subset$AI_prev)
    current_AI_prev_SD = sd(current_subset$AI_prev)
    current_Control_prev = mean(current_subset$Control_prev)
    current_Control_prev_SD = sd(current_subset$Control_prev)
    
    current_row = data.frame("Taxa_Name" = unique(current_subset$Taxa_Name),
                             "Iteration" = unique(current_subset$Iterations),
                             "LFC" = current_LFC, "LFC_SD" = current_LFC_SD,
                             "AI_Prev" = current_AI_prev, "AI_Prev_SD" = current_AI_prev_SD, 
                             "Control_Prev" = current_Control_prev, "Control_Prev_SD" = current_Control_prev_SD,
                             "Group" = current_method)
    
    current_data_Condensed = rbind(current_data_Condensed, current_row)
  }
  
  current_data = current_data_Condensed
  current_data$Taxa_Name = gsub("_", " ", current_data$Taxa_Name)
  long_data = current_data
  
  long_data$Trend = ifelse(long_data$LFC < 0, "Controls", 
                           ifelse(long_data$LFC > 0, "T1D", "NA"))
  
  Control_data = data.frame()
  T1D_data =data.frame()
  current_taxa = unique(long_data$Taxa_Name)[1]
  for (current_taxa in unique(long_data$Taxa_Name)) {
    current_taxa_subset = subset(long_data, long_data$Taxa_Name == current_taxa)
    if (current_taxa_subset$LFC > 0) {
      T1D_data = rbind(T1D_data, current_taxa_subset)
    } else {
      Control_data = rbind(Control_data, current_taxa_subset)
    }
  }
  
  if (current_method == "Reads/g") {
    current_method = "Reads/g"
  } else {
    current_method= "Relative Abundance"
  }
  
  
  LFC = ggplot(long_data, aes(x = LFC, y = Taxa_Name, color = Trend)) +
    geom_point( position = position_dodge(width=0.9)) +
    geom_errorbar(position = position_dodge(width=0.9),
                  aes(y=Taxa_Name, xmin=LFC-LFC_SD, xmax=LFC+LFC_SD),
                  width=0.4, alpha=0.9, size=1.3) +
    geom_vline(xintercept = 0) +
    ggtitle(paste(current_method, "ASV", sep = "\n")) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 7))

  assign(paste(current_method, current_glom, "plot", sep = "_"),LFC )
  print(paste(current_method, current_glom, "plot", sep = "_"))
}


jpeg("./Images/Deseq_Combined_All_ASVs.jpeg", res = 400, height = 4000, width = 5000)
ASV_plot = ggarrange(`Reads/g_Identifier_plot`, `Relative Abundance_Identifier_plot`, 
           ncol = 2,common.legend = T)

dev.off()
