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
setwd("~/Desktop/T1D_Demultiplexed/")
taxa_identifiers = read.csv("Taxa_identifiers_v2.csv")
Total_reads = readRDS("Filters_PS.RDS")
Sample_data = read.csv("T1D_Samples.csv")

#Matching list from Step 2 
Subset_1 = c("Region","Siblings_at_birth", "Apartment", "Total_Breastfeeding_Binned")

alpha = 0.05
Control_Num = 2
seed_number = 1
seed_max = 100

current_glom_1 = "Genus"
current_method = "Reads"
for (current_glom_1 in c("Genus", "Identifier")) {
  taxa_identifiers = read.csv("Taxa_identifiers_v2.csv")
  row.names(taxa_identifiers) = taxa_identifiers$ASV
  taxa_identifiers = taxa_identifiers[,c("ASV", current_glom_1)]
  
  if (current_glom_1 == "Genus") {
    ps.filtered.1 = tax_glom(Total_reads, "Genus")
    #Taxa Identification 
  } else {
    ps.filtered.1 = Total_reads
  }
  
  
  Actual_sign = data.frame()
  for (seed_number in 1:seed_max) {
    ps.filtered = ps.filtered.1
    if (seed_number%%10 == 0 ) {
      print(seed_number)
    }
    #####################
    #Sample Matching with n Controls
    subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
    subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")
    subset_controls = subset(subset_controls, subset_controls$ID != "20572")
    final_subset = subset_cases
    
    Matched = list()
    Count_1 = 0
    current_id = "2061"#unique(subset_df_1$ID)[1]
    for (current_id in unique(subset_cases$ID)){
      current_row = subset(subset_cases, subset_cases$ID == current_id)
      current_row_subset = current_row[,Subset_1]
      current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
      
      if (current_id == "2061") {
        current_row_subset = data.frame("Region"= current_row_subset)
      }
      
      Matched_subset = merge(subset_controls, current_row_subset)
      Matched = c(Matched, nrow(Matched_subset))
      #print(paste(current_id, nrow(Matched_subset)))
      
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
    #qPCR; 16s rRNA reads/g
    qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
    row.names(qpcr_data) = qpcr_data$ID
    qpcr_data$ID = NULL
    
    #Merge the Relative Abundance table and the qpcr table 
    otu_RA = data.frame(otu_table(ps.RA))
    otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
    row.names(otu_RA_1) = otu_RA_1$Row.names
    otu_RA_1$Row.names = NULL
    i = 4
    #Multiply all the Relabundance columns by the reads/g column
    for(i in 2:length(colnames(otu_RA_1))) {
      otu_RA_1[is.na(otu_RA_1[,i]),i] = 0
      otu_RA_1[,i] <- suppressWarnings(as.integer(round(otu_RA_1[,1] * otu_RA_1[, i])))
    }
    ps.RA.qpcr_1 = ps.RA
    otu_table(ps.RA.qpcr_1) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
    
    ps.current = ps.RA.qpcr_1
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
    count_table$ID = row.names(count_table)
    count_table_2 = merge(count_table, final_subset)
    AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
    Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
    count_table_2$Autoimmune_2_groups = factor(count_table_2$Autoimmune_2_groups, 
                                               levels = unique(count_table_2$Autoimmune_2_groups))
    
    OTU = data.frame(otu_table(ps.current))
    OTU[is.na(OTU)] = 0
    max(OTU)
    OTU <- mutate_all(OTU, function(x) as.integer(as.character(x)))
    otu_table(ps.current) = otu_table(OTU, taxa_are_rows = F)
    
    diagdds = suppressWarnings(phyloseq_to_deseq2(ps.current, ~ Autoimmune_2_groups))
    diagdds = estimateSizeFactors(diagdds, type = "poscounts")
    
    diagdds = DESeq(diagdds, test="Wald", fitType="local", quiet = TRUE)
    res = results(diagdds, cooksCutoff = FALSE)
    
    sigtab = data.frame(res[which(res$padj < alpha), ])
    
    if (nrow(sigtab)> 0) {
      sigtab = res[which(res$padj < alpha), ]
      sigtab = data.frame(sigtab)
      sigtab$ASV = rownames(sigtab)
      sigtab <- sigtab[,colSums(is.na(sigtab))<nrow(sigtab)]
      sigtab = merge(sigtab, taxa_identifiers)
      sigtab[,current_glom_1] = gsub("-", ".", sigtab[,current_glom_1])
      sigtab[,current_glom_1] = gsub("/", ".", sigtab[,current_glom_1])
      
      Identifier_list = unique(sigtab[,current_glom_1])

      current_taxa = Identifier_list[1]
      significant_data = data.frame()
      for (current_taxa in Identifier_list) {
        AI[,current_taxa] = as.numeric(AI[,current_taxa])
        Control[,current_taxa] = as.numeric(Control[,current_taxa])
        AI_prev = round(nrow(subset(AI, AI[,current_taxa] > 0))/nrow(AI), 3)*100
        Control_prev = round(nrow(subset(Control, Control[,current_taxa] > 0))/nrow(Control), 3)*100
        count_table_2[,current_taxa] = as.numeric(count_table_2[,current_taxa])
        count_table_2$Genera  = count_table_2[,current_taxa]
        P_value = subset(sigtab, sigtab[,current_glom_1] == current_taxa)
        P_value = P_value$padj
        current_row = data.frame("Taxa" = current_taxa, "AI_prev" = AI_prev, 
                                 "Control_prev" = Control_prev, "P_val" = P_value)
        significant_data = rbind(significant_data,current_row )
      }
    }
    sigtab = sigtab[,c(current_glom_1, "log2FoldChange", "baseMean")]
    colnames(sigtab) = c("Taxa_Name", "log2FoldChange", "baseMean" )
    Actually_significant = significant_data
    Actually_significant$Taxa_Name = Actually_significant$Taxa
    Actually_significant$Taxa = current_glom_1
    Actually_significant_1 = merge(Actually_significant, sigtab)
    Actual_sign = rbind(Actual_sign, Actually_significant_1)
  }
  Taxa_count = data.frame(table(Actual_sign$Taxa_Name))
  colnames(Taxa_count) = c("Taxa_Name", current_glom_1)
  assign(paste(current_glom_1, "_count", sep = ""), Taxa_count)
  write.csv(Taxa_count, paste("./CSV_files/Deseq_", current_method, "_", current_glom_1, "_Count.csv", sep = ""), row.names = F)
  write.csv(Actual_sign, paste("./CSV_files/Deseq_", current_method, "_", current_glom_1, "_2.csv", sep = ""), row.names = F)
}
