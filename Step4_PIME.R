#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr) #ggarrange
library(microbiomeMarker)

setwd("~/Desktop/T1D_Diabetologia/")
#####
# File import 
Sample_data = read.csv("Matched_Samples_Age.csv")
Sample_data$Group = ifelse(Sample_data$Autoimmune_2_groups == "Control", "Controls", "Type 1 Diabetes")
Sample_data$Age.collected.summary = paste("", Sample_data$Age.collected.summary, sep = "")
Sample_data_OG = Sample_data
rownames(Sample_data_OG) = Sample_data_OG$ID
subset_cases_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

rownames(Sample_data) = Sample_data$ID
Sample_data$Autoimmune_2_groups = ifelse(Sample_data$Autoimmune_2_groups == "Control", "Control", "Type 1 Diabetes")
Sample_data$Age.collected.summary = paste(Sample_data$Age.collected.summary, "_Month", sep = "")

####
#User input 
alpha = 0.05
seed_max = 10
Control_Num = 2

current_glom = "Identifier"
current_method = "Total"
for (current_method in c("Total", "Relative")) {
  for (current_glom in c("Genus", "Identifier")) {
    physeq_f = readRDS("Filters_PS.RDS") #Step 1
    physeq_tax = data.frame(tax_table(physeq_f))
    physeq_tax = subset(physeq_tax, is.na(physeq_tax$Genus) == F)
    tax_table(physeq_f) = tax_table(as.matrix(physeq_tax))
    sample_data(physeq_f) = sample_data(Sample_data_OG)
    
    if (current_glom == "Genus") {
      physeq_f = tax_glom(physeq_f, "Genus")
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
    
    ############################################
    #Begin
    ############################################

    PIME_MDA = data.frame()
    best_df = data.frame()
    OG_count = c()
    Count = c()
    for (seed_number in 1:seed_max) {
      if (current_method == "Total") {
        physeq_f = ps.RA.qpcr
      } else {
        physeq_f = ps.RA
      }
      
      ps_all = physeq_f
      set.seed(1)
      print(seed_number)
      subset_cases = subset_cases_OG
      subset_controls = subset_controls_OG
      PIME_data = subset_cases
      PIME_data$Group = "Type 1 Diabetes"
      
      set.seed(seed_number)
      Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num*16), ]
      Matched_subset$Group = "Controls"
      PIME_data = rbind(PIME_data,Matched_subset )
      
      table(PIME_data$Group)
      Sample_data = PIME_data
      Sample_data$Age.collected.summary = paste(Sample_data$Age.collected.summary, "", sep = "")
      rownames(Sample_data) = Sample_data$ID
      
      sample_data(physeq_f) = sample_data(Sample_data)
      physeq_f
      ps_all
      
      ################################################  
      count = c(count, pime.oob.error(physeq_f, "Group"))
      
      per_variable_obj= pime.split.by.variable(physeq_f, "Group")
      prevalences=pime.prevalence(per_variable_obj)
      set.seed(42)
      best.prev=pime.best.prevalence(prevalences, "Group")
      
      if (seed_number <= 10 ) {
        OG_count = c(OG_count, sum(phyloseq::sample_sums(physeq_f)))
        best.df = data.frame(best.prev$`OOB error`)
        best.df$Interval = NULL
        if (seed_number == 1) {
          best_df = best.df
        } else {
          if (nrow(best.df) == 18) {
            new_row = data.frame("OOB.error.rate...." = NA,
                                 "OTUs" = NA,"Nseqs" = NA)
            best.df = rbind(best.df, new_row)
          }
          best_df = cbind(best_df, best.df)
        }
      }
      
      
      if (current_glom == "Genus") {
        current_imp=best.prev$`Importance`$`Prevalence 70`
        prevalence.current= prevalences$`70`
        input_ord = ordinate(prevalence.current, "PCoA" , "binomial")
      } else {
        current_imp=best.prev$`Importance`$`Prevalence 50`
        prevalence.current= prevalences$`50`
        input_ord = ordinate(prevalence.current, "PCoA" , "binomial")
      }
      PIME_Tax = data.frame(tax_table(prevalence.current))
      PIME_Sample = data.frame(sample_data(prevalence.current))
      if (seed_number == 10) {
        graph_10 = plot_ordination(prevalence.current, input_ord , color = "Group", 
                                   shape = "Age.collected.summary")+
          stat_ellipse(aes(group = Group)) + 
          guides(shape=guide_legend(title="Month Stool\nCollected"),
                 color=guide_legend(title="Future\nDiagnosis")) +
          scale_shape_manual(breaks = c("6", "11", "12", "13", "14", "15", "17", "18", "NA"),
                             values = c(15,16,17, 23, 11,25, 7,4,8))+ 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          scale_color_manual(breaks = c("Controls", "Type 1 Diabetes"),
                             values = c("gray28", "dodgerblue")); graph_10
        
        #PIME no filtering 
        input_ord_0 = ordinate(physeq_f, "PCoA" , "binomial")
        graph_0 = plot_ordination(physeq_f, input_ord_0 , color = "Group", 
                                  shape = "Age.collected.summary")+
          stat_ellipse(aes(group = Group)) + 
          guides(shape=guide_legend(title="Month Stool\nCollected"),
                 color=guide_legend(title="Future\nDiagnosis")) +
          scale_shape_manual(breaks = c("6", "11", "12", "13", "14", "15", "17", "18", "NA"),
                             values = c(15,16,17, 23, 11,25, 7,4,8))+ 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          scale_color_manual(breaks = c("Controls", "Type 1 Diabetes"),
                             values = c("gray28", "dodgerblue")); graph_0
        #All no filtering 
        input_ord_all = ordinate(ps_all, "PCoA" , "binomial")
        graph_all = plot_ordination(ps_all, input_ord_all , color = "Group", 
                                    shape = "Age.collected.summary")+
          stat_ellipse(aes(group = Group)) + 
          guides(shape=guide_legend(title="Month Stool\nCollected"),
                 color=guide_legend(title="Future\nDiagnosis")) +
          scale_shape_manual(breaks = c("6", "11", "12", "13", "14", "15", "17", "18", "NA"),
                             values = c(15,16,17, 23, 11,25, 7,4,8))+ 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          scale_color_manual(breaks = c("Controls", "Type 1 Diabetes"),
                             values = c("gray28", "dodgerblue")); graph_all
        #ALL Filtered 
        
        if (current_glom == "Genus") {
          per_variable_obj_all= pime.split.by.variable(ps_all, "Group")
          prevalences_all=pime.prevalence(per_variable_obj_all)
          prevalence.current_all= prevalences_all$`70`
          input_ord_all= ordinate(prevalence.current_all, "PCoA" , "binomial")
        } else {
          per_variable_obj_all= pime.split.by.variable(ps_all, "Group")
          prevalences_all=pime.prevalence(per_variable_obj_all)
          prevalence.current_all= prevalences_all$`50`
          input_ord_all= ordinate(prevalence.current_all, "PCoA" , "binomial")
        }
        graph_all_filt = plot_ordination(prevalence.current_all, input_ord_all , color = "Group", 
                                         shape = "Age.collected.summary")+
          stat_ellipse(aes(group = Group)) + 
          guides(shape=guide_legend(title="Month Stool\nCollected"),
                 color=guide_legend(title="Future\nDiagnosis")) +
          scale_shape_manual(breaks = c("6", "11", "12", "13", "14", "15", "17", "18", "NA"),
                             values = c(15,16,17, 23, 11,25, 7,4,8))+ 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          scale_color_manual(breaks = c("Controls", "Type 1 Diabetes"),
                             values = c("gray28", "dodgerblue")); graph_all_filt
      }
      
      for (current_tax in unique(current_imp$SequenceID)) {
        current_subset = subset(current_imp, current_imp$SequenceID == current_tax)
        current_row = data.frame("ASV" = current_tax, 
                                 "MDA" = c(current_subset$Controls, 
                                           current_subset$Type.1.Diabetes),
                                 "Group" = c("Controls", "Type 1 Diabetes"))
        current_row = merge(current_row, taxa_identifiers)
        PIME_MDA = rbind(PIME_MDA, current_row)
      }#MDA loop 
    }#seed loop
    
    summary(count)
    PIME_Count = data.frame(table(PIME_MDA$Taxa))
    PIME_Count$Freq = PIME_Count$Freq/3
    PIME_Count = subset(PIME_Count, PIME_Count$Freq >= (seed_max/2))
    PIME_MDA_subset = subset(PIME_MDA, PIME_MDA$Taxa %in% PIME_Count$Var1)
    
    PIME_MDA_subset$Group = factor(PIME_MDA_subset$Group, 
                                   levels = c("Controls", "Type 1 Diabetes", "Total"))
    PIME_MDA_subset$Taxa = gsub("_NA_", "_", PIME_MDA_subset$Taxa )
    PIME_MDA_subset$Taxa = gsub("Clostridium_sensu_stricto_1", "Clostridium\nsensu stricto 1", PIME_MDA_subset$Taxa )
    
    for (i in 1:9) {
      PIME_MDA_subset$Taxa = gsub(paste("_", i, sep = ""), paste("\n", i, sep = ""), PIME_MDA_subset$Taxa )
    }
    PIME_MDA_subset$Taxa = gsub("_", " ", PIME_MDA_subset$Taxa )
    nrow(PIME_MDA_subset)
    current_tax = unique(PIME_MDA_subset$Taxa)[1]
    graph_data = data.frame()
    
    for (current_tax in unique(PIME_MDA_subset$Taxa)) {
      current_subset = subset(PIME_MDA_subset, PIME_MDA_subset$Taxa == current_tax)
      T1D = mean(subset(current_subset, current_subset$Group == "Type 1 Diabetes")$MDA)
      T1D_sd = sd(subset(current_subset, current_subset$Group == "Type 1 Diabetes")$MDA)
      
      Control = mean(subset(current_subset, current_subset$Group == "Controls")$MDA)
      Control_sd = sd(subset(current_subset, current_subset$Group == "Controls")$MDA)
      if (T1D < 0 | Control < 0) {
        PIME_MDA_subset = subset(PIME_MDA_subset, PIME_MDA_subset$Taxa != current_tax)
      } else {
        new_row = data.frame("Group" = c("Type 1 Diabetes", "Controls"),
                             "Average" = c(T1D, Control), 
                             "SD" = c(T1D_sd, Control_sd), 
                             "Taxa" = current_tax)
        graph_data = rbind(graph_data, new_row)
      }
    }
    final_plot = ggplot(PIME_MDA_subset, aes(x = Taxa, y = MDA, fill = Group)) +
      geom_hline(yintercept = 0) + geom_bar(position = "dodge", stat = "summary", fun = "mean")+
      ylab("Mean decrease accuracy (MDA)") + xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            legend.position="top",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_fill_manual(breaks = c("Controls", "Type 1 Diabetes" ),
                        values = c("#888888", "#6699CC"))
    
    graph_data$Taxa = gsub(" " , "\n", graph_data$Taxa)
    final_plot = ggplot(graph_data, aes(x = Taxa, y = Average, fill = Group)) +
      #geom_errorbar(aes(ymin = Average, ymax = Average+SD), position = position_dodge(0.9))+
      geom_bar(position = "dodge", stat = "identity")+
      ylab("Average(MDA)") + xlab("") + 
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            legend.position="top",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_fill_manual(breaks = c("Controls", "Type 1 Diabetes" ),
                        values = c("#888888", "#6699CC"))
    
    assign(paste(current_glom, "_MDA", sep = ""), final_plot)
    PCOA_final = ggarrange(graph_all, graph_all_filt, graph_0, graph_10, 
                           nrow = 4, common.legend = T, legend = "right")
    assign(paste(current_glom, "_PCOA", sep = ""), PCOA_final)
    
    jpeg(paste("PIME_MDA_", current_glom, "_", current_method, ".jpeg",sep = ""),
         res = 600, height = 2000, width = 4000)
    print(final_plot)
    dev.off()
    
    jpeg(paste("PIME_PCOA_", current_glom, "_", current_method,  ".jpeg", sep = ""), 
         res = 600, height = 4000, width = 2000)
    print(PCOA_final)
    dev.off()
    
    write.csv(graph_data, paste("PIME_iteration_", current_glom, "_", current_method,  ".csv", sep = ""), row.names = F)
  }#Genus or identifier
}
