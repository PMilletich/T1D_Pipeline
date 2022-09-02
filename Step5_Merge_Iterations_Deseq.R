library(ggplot2)
library(ggpubr)
setwd("~/Desktop/T1D_Demultiplexed/")

taxa_identifiers = read.csv("Taxa_Identifiers_v2.csv")
taxa_identifiers_1 = data.frame()
current_glom = "Identifier"
current_method = "Reads/g"
Final_Genus = c()
for (current_glom in c("Genus")) {
  iteration_Threshold = 50
  
  for (current_method in c("Reads/g", "Relative Abundance")) {
    if (current_method == "Reads/g") {
      current_data = read.csv(paste("./CSV_files/Deseq_Reads_", current_glom, "_2.csv", sep = ""))
    } else {
      current_data = read.csv(paste("./CSV_files/Deseq_RelAbun_", current_glom, "_2.csv", sep = ""))
    }
    
    df_Count = data.frame(table(current_data$Taxa_Name))
    df_Count = subset(df_Count, df_Count$Freq >= iteration_Threshold)
    
    Taxa_list = unique(df_Count$Var1)
    Taxa_list = subset(Taxa_list, grepl("NA.", Taxa_list) == F)
    
    Taxa_list_1 = c()
    if (current_glom == "Genus") {
      Final_Genus = append(Final_Genus, Taxa_list)
    } else {
      current_Genus = unique(Final_Genus)[1]
      for (current_Genus in unique(Final_Genus)) {
        current_list= subset(Taxa_list, grepl(current_Genus, Taxa_list) == T)
        if (current_Genus == "UCG.003") {
          current_list= subset(current_list, grepl("Erysipelotrichaceae_UCG.003", current_list) == F)
        }
        Taxa_list_1 = append(Taxa_list_1, current_list)
      }
      #Subset Genus list with final genus
      Taxa_list = Taxa_list_1
    }

    current_data = subset(current_data, current_data$Taxa_Name %in% Taxa_list)
    df_Count = data.frame(table(current_data$Taxa_Name))
    if (nrow(df_Count) > 0) {
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
    LFC = ggplot(long_data, aes(x = LFC, y = Taxa_Name, color = Trend)) +
      geom_point( position = position_dodge(width=0.9)) +
      geom_errorbar(position = position_dodge(width=0.9),
                    aes(y=Taxa_Name, xmin=LFC-LFC_SD, xmax=LFC+LFC_SD),
                    width=0.4, alpha=0.9, size=1.3) +
      geom_vline(xintercept = 0) +
      ggtitle(paste(current_method, current_glom, sep = "\n")) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(size = 7))



    # LFC_Controls = ggplot(Control_data, aes(x = LFC, y = Taxa_Name)) + 
    #   geom_point( position = position_dodge(width=0.9)) + 
    #   geom_errorbar(position = position_dodge(width=0.9),
    #                 aes(y=Taxa_Name, xmin=LFC-LFC_SD, xmax=LFC+LFC_SD), 
    #                 width=0.4, alpha=0.9, size=1.3) +
    #   geom_vline(xintercept = 0) +
    #   theme(axis.title.y = element_blank(), 
    #         axis.text.y = element_text(size = 7)) +
    #   ggtitle("Log2FoldChange; Greater in Controls")
    # 
    # 
    # if (nrow(T1D_data > 0)) {
    #   LFC_T1D = ggplot(T1D_data, aes(x = LFC, y = Taxa_Name)) + 
    #     geom_point( position = position_dodge(width=0.9)) + 
    #     geom_errorbar(position = position_dodge(width=0.9),
    #                   aes(y=Taxa_Name, xmin=LFC-LFC_SD, xmax=LFC+LFC_SD), 
    #                   width=0.4, alpha=0.9, size=1.3) +
    #     geom_vline(xintercept = 0) +
    #     theme(axis.title.y = element_blank(), 
    #           axis.text.y = element_text(size = 7)) +
    #     ggtitle("Log2FoldChange; Greater in T1D")
      
    #   LFC = ggarrange(LFC_Controls, LFC_T1D, ncol = 2)
    #   LFC = annotate_figure(LFC, paste(current_glom, current_method) )
    # } else {
    #   LFC = LFC_Controls
    #   LFC = annotate_figure(LFC, paste(current_glom, current_method) )
    # }

    if (current_method == "Reads/g") {
      current_method = "Reads"
    } else {
      current_method= "RelAbun"
    }
    jpeg(paste("./Images/Deseq_50_", current_method, "_", current_glom, ".jpeg", sep = ""), res = 400, height= 3000, width = 3000)
    print(LFC)
    dev.off()
    assign(paste(current_method, current_glom, "plot", sep = "_"),LFC )
    print(paste(current_method, current_glom, "plot", sep = "_"))
    }
  }
}

jpeg("./Images/Deseq_Combined_AG.jpeg", res = 400, height = 3000, width = 3000)
ggarrange(Reads_Genus_plot, RelAbun_Genus_plot,
          `Reads/g_Identifier_plot`, `Relative Abundance_Identifier_plot`,
          ncol = 2, nrow = 2, heights =c(1,1.5), common.legend = TRUE)
dev.off()



Greads_Count = read.csv("CSV_files/Deseq_Reads_Genus_Count.csv")
colnames(Greads_Count) = c("Taxa_Name", "Reads")
GRelAbun_Count = read.csv("CSV_files/Deseq_RelAbun_Genus_Count.csv")
colnames(GRelAbun_Count) = c("Taxa_Name", "RelAbun")

Genus_Count = merge(Greads_Count, GRelAbun_Count, all = T, by = "Taxa_Name")
Genus_Count = subset(Genus_Count, Genus_Count$Reads > 50 | Genus_Count$RelAbun > 50 )

ASVreads_Count = read.csv("CSV_files/Deseq_Reads_Identifier_Count.csv")
colnames(ASVreads_Count) = c("Taxa_Name", "Reads")
ASVRelAbun_Count = read.csv("CSV_files/Deseq_RelAbun_Identifier_Count.csv")
colnames(ASVRelAbun_Count) = c("Taxa_Name", "RelAbun")

ASV_Count = merge(ASVreads_Count, ASVRelAbun_Count, all = T)
ASV_Count = subset(ASV_Count, ASV_Count$Reads > 50 | ASV_Count$RelAbun > 50 )
ASV_Count= subset(ASV_Count, grepl("NA.", ASV_Count$Taxa_Name) == F)
