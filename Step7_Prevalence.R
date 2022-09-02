library(ggplot2)
library(ggforce)
library(vegan)
library(ggpubr)
library(DESeq2)
library(phyloseq)
library(dplyr)
'%notin%' = Negate('%in%')

###### 
#File import 
setwd("~/Desktop/T1D_V3/")
taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")
Total_reads = readRDS("Filters_PS.RDS")
Total_reads = tax_glom(Total_reads, "Genus")
Sample_data = read.csv("T1D_Samples.csv")

alpha = 0.05
Control_Num = 2
seed_number = 1
seed_max = 100


Prev_Difference = 25

current_glom_1 = "Genus"

current_method = "Matched"
Actual_sign = data.frame()
for (seed_number in 1:seed_max) {
  if (seed_number%%10 == 0 ) {
    print(seed_number)
  }
  ps.filtered = Total_reads
  #####################
  #Sample Matching with n Controls
  subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
  subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")
  
  final_subset = subset_cases
  Subset_1 = c("Region","Siblings_at_birth", "Apartment")
  
  Matched = list()
  Count_1 = 0
  current_id = "19289"#unique(subset_df_1$ID)[1]
  for (current_id in unique(subset_cases$ID)){
    current_row = subset(subset_cases, subset_cases$ID == current_id)
    current_row_subset = current_row[,Subset_1]
    current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
    
    if (current_id == "2061") {
      current_row_subset = data.frame("Region"= current_row_subset)
    }
    
    Matched_subset = merge(subset_controls, current_row_subset)
    Matched = c(Matched, nrow(Matched_subset))
    
    if (nrow(Matched_subset) >= 1 & typeof(current_row_subset) == "list") {
      set.seed(seed_number)
      Matched_subset = Matched_subset[sample(nrow(Matched_subset), Control_Num), ]
      final_subset = rbind(final_subset, Matched_subset)
      Count_1 = Count_1 + 1
    } else {
      set.seed(seed_number)
      Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num), ]
      final_subset = rbind(final_subset, Matched_subset)
    }
    subset_controls = subset(subset_controls, ! subset_controls$ID %in% final_subset$ID)
  }
  
  final_subset = unique(final_subset)
  table(final_subset$Autoimmune_Disorder)
  
  rownames(final_subset) = final_subset$ID
  sample_data(ps.filtered) = sample_data(final_subset)
  
  #################################################################
  #Relative Abundance; https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/transform_sample_counts
  ps.RA = transform_sample_counts(ps.filtered, function(x) x / sum(x) )
  #################################################################
  #qPCR 
  #################################################################
  qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
  rownames(qpcr_data) = qpcr_data$ID
  otu_RA = data.frame(otu_table(ps.RA))
  otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
  row.names(otu_RA_1) = otu_RA_1$Row.names
  otu_RA_1$Row.names = NULL
  otu_RA_1$ID = NULL
  otu_RA_2 = otu_RA_1
  
  for(i in 2:length(names(otu_RA_2))) {
    otu_RA_2[,i] <- round(otu_RA_2[,1] * otu_RA_2[, i])
  }
  otu_RA_2$copies_16s_per_gram_stool = NULL
  
  ps.RA.qpcr = ps.RA
  otu_RA_2[is.na(otu_RA_2)] = 0
  otu_RA_2 = otu_RA_2[rowSums(otu_RA_2[])>0,]
  otu_RA_2 = otu_RA_2[,colSums(otu_RA_2[])>0]
  
  otu_table(ps.RA.qpcr) = otu_table(otu_RA_2, taxa_are_rows = FALSE)
  ps.RA.qpcr
  
  #Taxa Identification 
  row.names(taxa_identifiers) = taxa_identifiers$ASV
  taxa_identifiers[,current_glom_1] = gsub("_V2", "", taxa_identifiers[,current_glom_1])
  taxa_identifiers[,current_glom_1] = gsub("_ASV", ".", taxa_identifiers[,current_glom_1])
  taxa_identifiers = taxa_identifiers[,c("ASV", current_glom_1)]
  
  ps.current = ps.RA.qpcr
  
  ###########
  #Data Fluffing
  count_table = data.frame(t(otu_table(ps.current)))
  count_table$ASV = row.names(count_table)
  count_table = merge(count_table, taxa_identifiers)
  row.names(count_table) = count_table[,current_glom_1]
  Identifier_list = unique(rownames(count_table))
  tax_table_all_group = unique(count_table[,current_glom_1])
  count_table$ASV = NULL
  count_table[,current_glom_1]= NULL
  count_table = data.frame(t(count_table))
  row.names(count_table) = gsub("X", "", row.names(count_table))
  Genus_list= colnames(count_table)
  count_table$ID = row.names(count_table)
  count_table_2 = merge(count_table, final_subset)
  AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
  Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
  count_table_2$Autoimmune_2_groups = factor(count_table_2$Autoimmune_2_groups, 
                                             levels = unique(count_table_2$Autoimmune_2_groups))
  current_genus = Genus_list[1]
  
  ttest_df = data.frame()
  current_genus = "Bifidobacterium"
  for (current_genus in Genus_list) {
    Trend = "No Trend"
    AI_prevalence = round(nrow(subset(AI, AI[,current_genus] > 0))/nrow(AI),3)*100
    Control_prevalence = round(nrow(subset(Control, Control[,current_genus] > 0))/nrow(Control),3)*100
    
    if (AI_prevalence == 0 & Control_prevalence > 0) {
      Trend = "Absent \nin T1D"
    } else if (AI_prevalence > 0 & Control_prevalence == 0 ) {
      Trend = "Absent \nin Control"
    } else if (abs(AI_prevalence - Control_prevalence) >= Prev_Difference ) {
      Trend = paste("Diff. \n>", Prev_Difference, sep = "")
    }
    
    if (Trend != "No Trend") {
      current_row = data.frame("Genus" = current_genus, 
                               "Prevalence" = AI_prevalence, 
                               "Group" = "T1D", 
                               "Trend" = Trend)
      Actual_sign = rbind(Actual_sign, current_row)
      current_row = data.frame("Genus" = current_genus, 
                               "Prevalence" = Control_prevalence, 
                               "Group" = "Control", 
                               "Trend" = Trend)
      Actual_sign = rbind(Actual_sign, current_row)
    } #Is there a Trend
  } #All Genus
} #Seed Number

Count = subset(Actual_sign, Actual_sign$Group == "Control")
Iterations = data.frame(table(Count$Genus))
Iterations = subset(Iterations, Iterations$Freq >= 50)

Actual_sign_Subset = subset(Actual_sign, Actual_sign$Genus %in% Iterations$Var1)
Actual_sign_Subset$Genus = gsub("_group", " group", Actual_sign_Subset$Genus)
Actual_sign_Subset$Genus = gsub("_", "\n", Actual_sign_Subset$Genus)
#Actual_sign_Subset = subset(Actual_sign_Subset, Actual_sign_Subset$Genus %in% Compare_Genus_list)

Actual_sign_Subset_Count = subset(Actual_sign_Subset, Actual_sign_Subset$Group == "Control")

Actual_sign_Subset_Count = Actual_sign_Subset_Count[order(Actual_sign_Subset_Count$Trend, Actual_sign_Subset_Count$Genus),]
Actual_sign_Subset_Count$Genus = factor(Actual_sign_Subset_Count$Genus, levels = unique(Actual_sign_Subset_Count$Genus))

Prev_plot = ggplot(Actual_sign_Subset, aes(x = Genus, y = Prevalence, color = Group)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6)) + 
  geom_boxplot() + 
  facet_grid(~Trend, scale = "free", space = "free")

jpeg("./Images/Prevalence_graph_Genus.jpeg", res = 400, height = 3000, width= 5000)
Prev_plot
dev.off()

unique(Actual_sign_Subset$Genus)
CD = subset(Actual_sign_Subset, Actual_sign_Subset$Genus == "Candidatus\nSoleaferrea")
CD = subset(CD, CD$Group != "T1D")
mean(CD$Prevalence)
