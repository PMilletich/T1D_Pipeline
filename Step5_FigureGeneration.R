#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr) #ggarrange
library(microbiomeMarker)

setwd("~/Desktop/T1D_Diabetologia/")
#####
# File import 
physeq_f = readRDS("Filters_PS.RDS") #Step 1
physeq_tax = data.frame(tax_table(physeq_f))
physeq_tax = subset(physeq_tax, is.na(physeq_tax$Genus) == F)
tax_table(physeq_f) = tax_table(as.matrix(physeq_tax))

Sample_data = read.csv("Matched_Samples.csv")
Sample_data$Group = ifelse(Sample_data$Autoimmune_2_groups == "Control", "Controls", "Type 1 diabetes")
Sample_data$Age.collected.summary = paste("", Sample_data$Age.collected.summary, sep = "")
Sample_data_OG = Sample_data
rownames(Sample_data_OG) = Sample_data_OG$ID

subset_cases_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

rownames(Sample_data) = Sample_data$ID
Sample_data$Autoimmune_2_groups = ifelse(Sample_data$Autoimmune_2_groups == "Control", "Control", "Type 1 diabetes")
Sample_data$Age.collected.summary = paste(Sample_data$Age.collected.summary, "_Month", sep = "")

sample_data(physeq_f) = sample_data(Sample_data_OG)

####
#User input 
alpha = 0.05
seed_max = 100
Control_Num = 2
prev_thresh = 0.1
current_method = "Relative"
current_glom = "Genus"

#######
#Relative and Total Abundance 
#################################################################
#Abundance 
#################################################################
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

qpcr_data = Sample_data_OG[,c("ID", "copies_16s_per_gram_stool")]
rownames(qpcr_data) = qpcr_data$ID
qpcr_data$ID = NULL

otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL

for(i in 2:length(names(otu_RA_1))) {
  otu_RA_1[,i] <- round(otu_RA_1[,1] * otu_RA_1[, i])
}

otu_RA_1$copies_16s_per_gram_stool = NULL

ps.RA.qpcr = ps.RA
otu_table(ps.RA.qpcr) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
ps.RA.qpcr

for (current_method in c("Total", "Relative")) {
  for (current_glom in c("Genus", "Identifier")) { #
    if (current_method == "Total") {
      ps_current_1 = ps.RA.qpcr
    } else {
      ps_current_1 = ps.RA
    }
    
    
    #######
    #Deseq
    
    Actual_sign_subset = read.csv(paste("./CSV_files/Deseq_", current_method, "_", 
                                 current_glom, "LFC.csv",sep = ""))
    Actual_sign_subset$Group = gsub("Diabetes", "diabetes", Actual_sign_subset$Group)
    
    foo <- function(x) ceiling(x/10)*10
  
    max = 30
    min = -30

    LFC_plot = ggplot(Actual_sign_subset, aes(x = Taxa, y = log2FoldChange, color = Group)) + 
      geom_boxplot() + xlab("") + 
      scale_y_continuous(limits = c(min, max),
                         breaks = c(min, median(seq(min, 0)), 0, median(seq(0,max)), max)) + 
      geom_hline(yintercept = 0, alpha = 0.5) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            legend.position="top",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      scale_color_manual(breaks = c("Controls", "Type 1 diabetes" ),
                         values = c("#888888", "#6699CC")); LFC_plot
    
    assign(paste(current_glom, current_method, "Deseq", sep = "_"), LFC_plot)
    
    
    if (current_glom == "Genus") {
      ps_current_1 = tax_glom(ps_current_1, "Genus")
      
      taxa_identifiers = read.csv("Taxa_identifiers_v2.csv")
      row.names(taxa_identifiers) = taxa_identifiers$ASV
      taxa_identifiers = taxa_identifiers[,c("ASV", "Genus")]
      colnames(taxa_identifiers) = c("ASV", "Taxa")
    } else {
      taxa_identifiers = read.csv("Taxa_identifiers_v2.csv")
      row.names(taxa_identifiers) = taxa_identifiers$ASV
      taxa_identifiers = taxa_identifiers[,c("ASV", "Identifier")]
      colnames(taxa_identifiers) = c("ASV", "Taxa")
    }
    
    Taxa_list = read.csv(paste("./Review_2/PIME_iteration_", current_glom, 
                               "_", current_method, ".csv", sep = ""))
    Taxa_list = sort(unique(Taxa_list$Taxa))
    MDA_data = read.csv(paste("./Review_2/PIME_iteration_", current_glom, 
                              "_", current_method, ".csv", sep = ""))
    
    ############################################
    #Begin
    ############################################
    sign_data = data.frame()
    for (seed_number in 1:seed_max) {
      ps_current= ps_current_1
      subset_cases = subset_cases_OG
      subset_controls = subset_controls_OG
      final_data = subset_cases
      final_data$Group = "Type 1 diabetes"
      
      set.seed(seed_number)
      Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num*16), ]
      Matched_subset$Group = "Controls"
      final_data = rbind(final_data,Matched_subset )
      
      table(final_data$Group)
      final_data$Age.collected.summary = paste(final_data$Age.collected.summary, "", sep = "")
      rownames(final_data) = final_data$ID
      
      sample_data(ps_current) = sample_data(final_data)
      
      ################################################  
      #Data fluff
      count_table = data.frame(t(otu_table(ps_current)))
      count_table$ASV = row.names(count_table)
      count_table = merge(count_table, taxa_identifiers)
      row.names(count_table) = count_table$Taxa
      Identifier_list = unique(rownames(count_table))
      tax_table_all_group = unique(count_table$Taxa)
      count_table$ASV = NULL
      count_table$Taxa= NULL
      count_table = data.frame(t(count_table))
      row.names(count_table) = gsub("X", "", row.names(count_table))
      count_table$ID = row.names(count_table)
      count_table_2 = merge(count_table, final_data)
      AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
      Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
      count_table_2$Autoimmune_2_groups = factor(count_table_2$Autoimmune_2_groups, 
                                                 levels = unique(count_table_2$Autoimmune_2_groups))
      
      Average_Data = data.frame()
      current_taxa = Taxa_list[1]
      for (current_taxa in Taxa_list) {
        if (current_glom == "Identifier") {
          current_taxa_split = strsplit(current_taxa, "\n") 
          Genus = current_taxa_split[[1]][1]
          ID_num = paste("_", tail(current_taxa_split[[1]][length(current_taxa_split[[1]])]), sep = "")
          
          col_X = subset(colnames(AI), grepl(ID_num, colnames(AI)) == T)
          col_X = subset(col_X, grepl(Genus, col_X) == T)
          AI_count = AI[,col_X]
          Control_count = Control[,col_X]
        } else {
          AI_count = AI[,current_taxa]
          Control_count = Control[,current_taxa]
        }

        AI_mean = mean(AI_count)
        Control_mean = mean(Control_count)
        
        Difference = AI_mean - Control_mean
        
        New_data = data.frame("Difference" = Difference,
                              "Taxa" = current_taxa)
        Average_Data= rbind(Average_Data, New_data)
      } #Taxa loop 
      sign_data = rbind(sign_data, Average_Data)
    } #Seed loop 
    
    sign_data$Group = ""
    for (current_tax in unique(sign_data$Taxa)) {
      current_subset = subset(sign_data, sign_data$Taxa == current_tax)
      if (mean(current_subset$Difference) > 0) {
        sign_data[sign_data$Taxa == current_tax,]$Group = "Type 1 diabetes"
      } else {
        sign_data[sign_data$Taxa == current_tax,]$Group = "Controls"
      }
    }
    
    foo_1 <- function(x) ceiling(max(x)/1000)*1000
    foo_1M <- function(x) (ceiling(x/10)-10)*10
    if (current_method == "Total") {
      
      max = foo_1(max(sign_data$Difference))
      min = foo_1M(min(sign_data$Difference))
      break_range = c(round(min,-3), round(median(seq(min, 0)),-3), 0, round(median(seq(0, max)),-3), round(max,-3))
    } else {
      
      max = max(sign_data$Difference)
      min = min(sign_data$Difference)
      break_range = c(round(min,3), round(median(seq(min, 0, 0.001)),3), 0, round(median(seq(0, max, 0.001)),3), round(max,3))
    }
    
    
    Mean_Abundance = ggplot(sign_data, aes(x = Taxa, y = Difference, color = Group)) + 
      geom_boxplot() + 
      xlab("") + 
      scale_y_continuous(name="Difference in mean abundance", labels = scales::comma,
                         limits = c(min, max),
                         breaks = break_range) + 
      #scale_y_discrete(position = "right") +
      geom_hline(yintercept = 0, alpha = 0.25) +
      theme(axis.text.x = element_text(angle = 45, hjust=1),
        #axis.text.y=element_blank(),    
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_color_manual(breaks = c("Controls", "Type 1 diabetes"),
                         values = c("#888888", "#6699CC")); Mean_Abundance
   
    
    MDA_data$Group = gsub("Diabetes", "diabetes", MDA_data$Group)
    dodge = position_dodge(width=0.9)
    
    max = max(MDA_data$Average)
    
    Mean_MDA = ggplot(MDA_data, aes(x = Taxa, y = Average, fill = Group)) +
       geom_bar(stat = "identity", position = dodge) +
       scale_y_continuous(name="Average(MDA)",
                          limits = c(0, max),
                          breaks = c(0, round(median(seq(0,max, 0.001)),3), round(max,2))) + 
       ylab("") + xlab("")+
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            #axis.text.y=element_blank(),
            #legend.position="none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_fill_manual(breaks = c("Controls", "Type 1 diabetes"),
                         values = c("#888888", "#6699CC")); Mean_MDA
     
     assign(paste(current_glom, "_", current_method, "_MDA_0", sep = ""), Mean_MDA)
     assign(paste(current_glom, "_", current_method, "_Abundance", sep = ""), Mean_Abundance)
     
     print(paste(current_glom, "_", current_method, sep = ""))
  }#Genus or Identifier 
}#Total or Relative 

jpeg('Fig2.jpeg', res = 650, height = 6000, width = 6000)
ggarrange(Genus_Total_MDA_0, Genus_Relative_MDA_0,
          Genus_Total_Abundance, Genus_Relative_Abundance,
          Genus_Total_Deseq, Genus_Relative_Deseq,
          nrow = 3, ncol = 2,common.legend = T, align = "hv")
dev.off()

jpeg('Fig3.jpeg', res = 650, height = 6000, width = 9000)
ggarrange(Identifier_Total_MDA_0, Identifier_Relative_MDA_0,
          Identifier_Total_Abundance, Identifier_Relative_Abundance,
          Identifier_Total_Deseq, Identifier_Relative_Deseq,
          nrow = 3, ncol = 2,common.legend = T, align = "hv")
dev.off()





