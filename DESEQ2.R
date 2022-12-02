library(ggplot2)
library(DESeq2)
library(dplyr) 

###### 
#File import 
taxa_identifiers = read.csv("Taxa_identifiers_v2.csv")
Total_reads = readRDS("Filters_PS.RDS") #Step 1
physeq_tax = data.frame(tax_table(Total_reads))
physeq_tax = subset(physeq_tax, is.na(physeq_tax$Genus) == F)
tax_table(Total_reads) = tax_table(as.matrix(physeq_tax))

Sample_data = read.csv("Matched_Samples.csv")
rownames(Sample_data) = Sample_data$ID
sample_data(Total_reads) = sample_data(Sample_data)

qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
row.names(qpcr_data) = qpcr_data$ID
qpcr_data$ID = NULL

#User Input 
alpha = 0.05
Control_Num = 2
seed_max = 100

#Starting values (Will be overwritten in loop)
seed_number = 1
current_method = "Relative"
current_glom_1 = "Genus"

#Convert Relative abundance to ints by multiplying smallest value 
ps.RA.ALL = transform_sample_counts(Total_reads, function(x) x / sum(x) )
OTU_All = data.frame(otu_table(ps.RA.ALL))
list = c()
for (i in 1:ncol(OTU_All)) {
  current_list = OTU_All[,i]
  list = c(list, current_list)
}
list= subset(list, list != 0)
summary(list)
OTU_All = OTU_All*1e6
otu_table(ps.RA.ALL) = otu_table(OTU_All, taxa_are_rows = F)

#Total Abundance 
#################################################################
#Relative Abundance; https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/transform_sample_counts
ps.RA = transform_sample_counts(Total_reads, function(x) x / sum(x) )

otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL
i = 4
for(i in 2:length(colnames(otu_RA_1))) {
  otu_RA_1[is.na(otu_RA_1[,i]),i] = 0
  otu_RA_1[,i] <- suppressWarnings(as.integer(round(otu_RA_1[,1] * otu_RA_1[, i])))
}
ps.RA.qpcr_1 = ps.RA
otu_table(ps.RA.qpcr_1) = otu_table(otu_RA_1, taxa_are_rows = FALSE)

for (current_method in c("Total", "Relative")) {
  for (current_glom_1 in c("Genus", "Identifier")) {
    if (current_method == "Total") {
      ps.filtered.1 = ps.RA.qpcr_1
    } else {
      ps.filtered.1 = ps.RA.ALL
    }
    
    taxa_identifiers = read.csv("Taxa_identifiers_v2.csv")
    row.names(taxa_identifiers) = taxa_identifiers$ASV
    taxa_identifiers = taxa_identifiers[,c("ASV", current_glom_1)]
    
    if (current_glom_1 == "Genus") {
      ps.filtered.1 = tax_glom(ps.filtered.1, "Genus")
      prev_thresh = 0.1
      #Taxa Identification 
    } else {
      ps.filtered.1 = ps.filtered.1
      prev_thresh = 0.2
    }
    
    Actual_sign = data.frame()
    Actual_sign_subset= data.frame()
    for (seed_number in 1:seed_max) {
      if (seed_number%%10 == 0) {
        print(seed_number)
      }
      ps.filtered = ps.filtered.1
      
      subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
      subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")
      subset_controls = subset(subset_controls, subset_controls$ID != "20572")
      final_subset = subset_cases
      
      set.seed(seed_number)
      Matched_subset = subset_controls[sample(nrow(subset_controls), (Control_Num*16)), ]
      final_subset = rbind(final_subset, Matched_subset)
      table(final_subset$Autoimmune_Disorder)
      rownames(final_subset) = final_subset$ID
      sample_data(ps.filtered) = sample_data(final_subset)
      
      ps.current = ps.filtered
      ###########
      #Data fluff
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
      #Prevalence filter
      ASV_list = colnames(OTU)
      OTU$ID = rownames(OTU)
      for (current_tax in ASV_list) {
        current_subset = OTU[,c("ID", current_tax)]
        T1D = subset(current_subset, current_subset$ID %in% subset_cases$ID)
        T1D_prev = round(nrow(subset(T1D, T1D[,current_tax] > 0))/nrow(T1D),3)
        Controls = subset(current_subset, ! current_subset$ID %in% subset_cases$ID)
        Controls_prev = round(nrow(subset(Controls, Controls[,current_tax] > 0))/nrow(Controls),3)
        if (T1D_prev < prev_thresh & Controls_prev < prev_thresh) {
          OTU[,current_tax] = NULL
        }
      }
      count_table$ID = NULL
      
      OTU[is.na(OTU)] = 0
      OTU <- mutate_all(OTU, function(x) as.integer(as.character(x)))
      otu_table(ps.current) = otu_table(OTU, taxa_are_rows = F)
      
      #DESEQ
      diagdds = suppressWarnings(phyloseq_to_deseq2(ps.current, ~ Autoimmune_2_groups))
      diagdds = estimateSizeFactors(diagdds, type = "poscounts")
      diagdds = DESeq(diagdds, test="Wald", fitType="local", quiet = TRUE)
      res = results(diagdds, cooksCutoff = FALSE)
      
      sigtab = data.frame(res[which(res$padj < alpha), ])
      sigtab$ASV = rownames(sigtab)
      sigtab= merge(sigtab, taxa_identifiers)
      
      if (nrow(sigtab)> 0) {
        sigtab[,current_glom_1] = gsub("-", ".", sigtab[,current_glom_1])
        Identifier_list = unique(sigtab[,current_glom_1])
        significant_data = sigtab[,c("log2FoldChange", current_glom_1)]
        Actual_sign = rbind(Actual_sign, significant_data)
      }
    }#Seed numbers 
    Taxa_count = data.frame(table(Actual_sign[,current_glom_1]))
    Taxa_count = subset(Taxa_count, Taxa_count$Freq >= (seed_max/2))
    
    Actual_sign_subset = subset(Actual_sign, Actual_sign[,current_glom_1] %in% Taxa_count$Var1)
    colnames(Actual_sign_subset) = c("log2FoldChange", "Taxa")
    Actual_sign_subset$Group = ""
    
    for (current_tax in unique(Actual_sign_subset$Taxa)) {
      current_subset = subset(Actual_sign_subset, Actual_sign_subset$Taxa == current_tax)
      Average = mean(current_subset$log2FoldChange)
      if (Average > 0) {
        current_Group = "Type 1 Diabetes"
      } else {
        current_Group = "Controls"
      }
      Actual_sign_subset[Actual_sign_subset$Taxa == current_tax,]$Group = current_Group
    }
    
    Actual_sign_subset$Taxa = gsub("_NA_", "_", Actual_sign_subset$Taxa)
    Actual_sign_subset$Taxa = gsub("Clostridium_sensu_stricto_1_", "Clostridium\nsensu stricto 1\n", Actual_sign_subset$Taxa)
    
    Actual_sign_subset$Taxa = gsub("_", "\n", Actual_sign_subset$Taxa)
    
    LFC_plot = ggplot(Actual_sign_subset, aes(x = Taxa, y = log2FoldChange, color = Group)) + 
      geom_boxplot() + xlab("") + 
      geom_hline(yintercept = 0, alpha = 0.5) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            legend.position="top",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      scale_color_manual(breaks = c("Controls", "Type 1 Diabetes" ),
                         values = c("#888888", "#6699CC")); LFC_plot
    
    jpeg(paste("Deseq_", current_glom_1, "_", prev_thresh, "_", current_method, "_TA.jpeg", sep = ""), 
         res = 600, height = 2000, width = 4000)
    print(LFC_plot)
    dev.off()
  }
}
