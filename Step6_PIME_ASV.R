#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/T1D_Demultiplexed/")

###############
#File Import 
physeq_f = readRDS("Filters_PS.RDS")
Sample_data = data.frame(sample_data(physeq_f))

taxa_identifiers = read.csv("Taxa_Identifiers_v2.csv")
row.names(taxa_identifiers) = taxa_identifiers$ASV
taxa_identifiers = taxa_identifiers[,c("ASV", "Identifier")]

###############cdxffffffffffffffffffffffffffffffffffffffffffffffffffffffff
#Default Values 

Control_Num = 2
Subset_1 = c("Region","Siblings_at_birth", "Apartment", "Total_Breastfeeding_Binned")
seed_max = 100
seed_number = 1
current_method = "Reads"

subset_cases_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

#################################################################
#qPCR 
#################################################################
qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
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

#################################################################
#Begin Running 
#################################################################
for (current_method in c("Reads", "RelativeAbundance")) {
  total_data = data.frame()
  count_data = data.frame()
  best_df = data.frame()
  for (seed_number in 1:seed_max) {
    print(seed_number)
    subset_cases = subset_cases_OG
    subset_controls = subset_controls_OG
    PIME_data = subset_cases
    PIME_data$Group = "T1D"
    
    current_id = "19289"#unique(subset_df_1$ID)[1]
    for (current_id in unique(subset_cases$ID)){
      current_row = subset(subset_cases, subset_cases$ID == current_id)
      current_row_subset = current_row[,Subset_1]
      current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
      
      if (current_id == "2061") {
        current_row_subset = data.frame("Region"= current_row_subset)
      }
      
      Matched_subset = merge(subset_controls, current_row_subset)
      if (nrow(Matched_subset) >= 1 & typeof(current_row_subset) == "list") {
        set.seed(seed_number)
        Matched_subset = Matched_subset[sample(nrow(Matched_subset), Control_Num), ]
        Matched_subset$Group = "Controls"
        PIME_data = rbind(PIME_data, Matched_subset)
        subset_controls = subset(subset_controls, ! subset_controls$ID %in% PIME_data$ID)
      } else {
        set.seed(seed_number)
        Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num), ]
        Matched_subset$Group = "Controls"
        PIME_data = rbind(PIME_data, Matched_subset)
        subset_controls = subset(subset_controls, ! subset_controls$ID %in% PIME_data$ID)
      }
    }
    table(PIME_data$Group)
    Sample_data = PIME_data
    rownames(Sample_data) = Sample_data$ID
    
    if (current_method == "Reads") {
      physeq_f = ps.RA.qpcr
    } else {
      physeq_f = ps.RA
    }
    sample_data(physeq_f) = sample_data(Sample_data)
    physeq_f
    ################################################
    #Celiac
    print(pime.oob.error(physeq_f, "Group"))
    per_variable_obj= pime.split.by.variable(physeq_f, "Group")
    
    prevalences=pime.prevalence(per_variable_obj)
    set.seed(42)
    best.prev=pime.best.prevalence(prevalences, "Group")
    
    if (seed_number <= 10 ) {
      print(paste(seed_number, sum(phyloseq::sample_sums(physeq_f))))
      best.df = data.frame(best.prev$`OOB error`)
      best.df$Iteration = seed_number
      if (seed_number == 1) {
        best_df = best.df
      } else {
        if (nrow(best.df) == 18) {
          new_row = data.frame("Interval" = "Prevalence 95%","OOB.error.rate...." = NA,
                               "OTUs" = NA,"Nseqs" = NA, "Iteration" = NA)
          best.df = rbind(best.df, new_row)
        }
        best_df = cbind(best_df, best.df)
      }
    }
    
    imp50=best.prev$`Importance`$`Prevalence 50`
    prevalence.50 = prevalences$`50`
    input_ord = ordinate(prevalence.50, "PCoA" , "binomial")
    
    if (seed_number == 10) {
      graph_10 = plot_ordination(prevalence.50, input_ord , color = "Group")+
        stat_ellipse(aes(group = Group)) + 
        ggtitle(paste("ASV\n", current_method, sep = ""))
    }
    colnames(taxa_identifiers) = c("SequenceID", "Identifier")
    imp50_subset = subset(imp50, imp50$MeanDecreaseAccuracy > 0 &
                            imp50$Controls > 0 &
                            imp50$T1D > 0 )
    imp50_subset = merge(imp50_subset, taxa_identifiers)
    ###########
    #Data Fluffing
    colnames(taxa_identifiers) = c("ASV", "Identifier")
    count_table = data.frame(t(otu_table(prevalence.50)))
    count_table$ASV = row.names(count_table)
    count_table = merge(count_table, taxa_identifiers)
    row.names(count_table) = count_table$Identifier
    Identifier_list = unique(rownames(count_table))
    tax_table_all_group = unique(count_table$Identifier)
    count_table$ASV = NULL
    count_table$Identifier= NULL
    count_table = data.frame(t(count_table))
    row.names(count_table) = gsub("X", "", row.names(count_table))
    count_table$ID = row.names(count_table)
    count_table_2 = merge(count_table, PIME_data)
    AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
    Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
    
    current_Identifier = unique(imp50_subset$Identifier)[1]
    mean_data = data.frame()
    for (current_Identifier in unique(imp50_subset$Identifier)) {
      current_subset = subset(imp50_subset, imp50_subset$Identifier == current_Identifier)
      current_subset = current_subset[,c("Controls","T1D","MeanDecreaseAccuracy","Identifier")]
      count_data = rbind(count_data, current_subset)
      
      current_Identifier = gsub("/", ".", current_Identifier)
      current_Identifier = gsub("-", ".", current_Identifier)
      
      AI_mean_df = data.frame("Values" = as.numeric(AI[,current_Identifier]),
                              "Genus" = current_Identifier, 
                              "Group" = "T1D", 
                              "MDA" = current_subset$T1D)
      
      
      Control_mean_df = data.frame("Values" = as.numeric(Control[,current_Identifier]),
                                   "Genus" = current_Identifier, 
                                   "Group" = "Control", 
                                   "MDA" = current_subset$Controls)
      
      current_subset = rbind(AI_mean_df, Control_mean_df)
      mean_data= rbind(mean_data, current_subset)
    }
    total_data = rbind(total_data, mean_data)
  }
  assign(paste(current_method, "total", sep = "_"), total_data)
  assign(paste(current_method, "count", sep = "_"), count_data)
  assign(paste("PCOA", current_method, sep = "_"), graph_10)
}


final_graph = ggarrange(PCOA_Reads, PCOA_RelativeAbundance, ncol = 2, 
                        common.legend = T)
jpeg(paste("./Images/PCOA_PIME_50_ASV.jpeg", sep = ""), res = 400, height = 2000, width = 4000)
final_graph
dev.off()

#write.csv(RelativeAbundancetotal, "./CSV_Output/RelAbund_PIME_100.csv", row.names = F)
#write.csv(Readstotal, "./CSV_Output/Reads_PIME_100.csv", row.names = F)

current_method = "Reads/g"
for (current_method in c("Reads/g", "Relative Abundance")) {
  if (current_method == "Reads/g") {
    total_data = Reads_total
    count_data = Reads_count
    mean_num = 1
  } else {
    total_data = RelativeAbundance_total 
    count_data = RelativeAbundance_count
    mean_num = 5
  }
  
  count_data = data.frame(table(count_data$Identifier))
  count_data = subset(count_data, count_data$Freq >= 50)
  
  total_data_subset = subset(total_data, total_data$Genus %in% count_data$Var1)
  total_data_subset = subset(total_data_subset, grepl("NA.", total_data_subset[,2]) == F)
  
  total_data_T1D = subset(total_data_subset, total_data_subset$Group == "T1D")
  total_data_Control = subset(total_data_subset, total_data_subset$Group != "T1D")
  
  T1D_data = data.frame()
  Control_data = data.frame()
  current_genus = unique(total_data_subset$Genus)[1]
  for (current_genus in unique(total_data_subset$Genus)) {
    AI_subset = subset(total_data_T1D, total_data_T1D$Genus == current_genus)
    AI_mean = mean(AI_subset[,1])
    
    Control_subset = subset(total_data_Control, total_data_Control$Genus == current_genus)
    Control_mean = mean(Control_subset[,1])
    if (AI_mean > Control_mean) {
      T1D_data = rbind(T1D_data,AI_subset)
      T1D_data = rbind(T1D_data,Control_subset)
    } else if (AI_mean == Control_mean) {
      print(paste("Uh oh", current_genus))
    } else {
      Control_data = rbind(Control_data,AI_subset)
      Control_data = rbind(Control_data,Control_subset)
    }
  }
  
  T1D_data$Genus = gsub("_", "\n", T1D_data$Genus)
  Control_data$Genus = gsub("_", "\n", Control_data$Genus)
  T1D_graph = ggplot(T1D_data,  aes(x = Genus, y = Values, color = Group)) + 
    geom_boxplot(outlier.alpha = 0.5) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank()) + ylab(current_method) +
    ggtitle(paste(current_method, "\nGreater in T1D")) +
    coord_cartesian(ylim = c(0, mean(c(T1D_data$Values, Control_data$Values), na.rm = T)*mean_num))
  
  Control_graph = ggplot(Control_data,  aes(x = Genus, y = Values, color = Group)) + 
    geom_boxplot( outlier.alpha = 0.5) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank()) + ylab(current_method) + 
    ggtitle(paste(current_method, "\nGreater in Controls"))  +
    coord_cartesian(ylim = c(0, mean(c(T1D_data$Values, Control_data$Values), na.rm = T)*mean_num))
  
  box_graph = ggarrange(T1D_graph, Control_graph, nrow = 2, common.legend = T)
  box_graph
  assign(paste(current_method, "graph", sep = "_"), box_graph)
}
plot_final = ggarrange(`Reads/g_graph`,`Relative Abundance_graph`,
          ncol =2 , common.legend = T); plot_final

jpeg ("./Images/PIME_RR_ASV_50.jpeg", res = 400, height = 4000, width = 4000)
plot_final
dev.off()

