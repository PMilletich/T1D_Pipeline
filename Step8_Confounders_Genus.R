
#####################
#Libraries and Functions 
#####################
library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
library(ggplot2)
library(ggpubr)

'%notin%' = Negate('%in%')#https://stackoverflow.com/questions/38351820/negation-of-in-in-r
setwd("~/Desktop/T1D_Demultiplexed/")
#####################
#Load Files 
#####################
taxa_identifiers = read.csv("Taxa_Identifiers_v2.csv")

#####################
#Default Values 
#####################
#Create list of columns to compare 

environmental_list = c("DR13.DQ603","DR13.DQ604","DR15.DQ602","DR3.DQ2.5",
                       "DR4.DQ7","DR4.DQ8","DR5.DQ7","DR7.DQ2","DR7.DQ9","DR8.DQ4","DR9.DQ9","DR1.DQ5",
                       "Antibiotics_during_Pregnancy","Infection_during_Pregnancy","Infection_as_newborn","Infection_w_antibiotics_1.12m",
                       "Cold.UpperRespTractInf","Otitis","Pneumonia","gastroenteritis","OtherInfection","OtherDisease",
                       "Corticosteroids_Pregnancy","HighBloodPressure_Meds_Pregnancy",
                       "PainKiller_Meds_Pregnancy","HormonePreparates_Pregnancy","Other_Meds_Pregnancy","Total_Meds_Pregnancy",
                       "Intro_CowsMilk_Binned","Intro_Gluten_Binned","Egg_first_year","Beef_first_year","Pork_first_year",
                       "Exclusive_Breastfeeding_Binned","Total_Breastfeeding_Binned","Intro_formula_Binned",
                       #"Autoimmune_2_groups","Region","Sex","Siblings_at_birth",
                       "Delivery","Apartment","Father_Only_Elementary",
                       "Father_Unemployed","Mother_Unemployed","Both_Parents_Abroad","Single_Mother",
                       "Stressful_Life_Event",
                       "Worry_ChronicIllness_Child", "Smoking_Pregnancy","Risky_Alcohol_Pregnancy",
                       #"Alcohol_Pregnancy_PerMonth",
                       "Alcohol_Smoking_Medications_Pregnancy","Mother_Over_35","Father_Over_40" )

#Create lists of taxa to analyze
Deseq_Genus_Control = c("Parasutterella")
Deseq_Genus_T1D = c()

PIME_Genus_T1D = c("Alistipes", "Bacteroides", "Fusicatenibacter", 
                   "Hungatella", "Roseburia","Ruminococcus")
PIME_Genus_Control = c("Dorea", "Flavonifractor", "UBA1819")


table(map$Risky_Alcohol_Pregnancy)
for (current_glom in c("Genus")) {
  PS = readRDS("Filters_PS.RDS")
  map = read.csv("Matched_Samples.csv")
  rownames(map) = map$ID
  sample_data(PS) = sample_data(map)
  
  if (current_glom == "Genus") {
    #Create lists of Taxa of interests from steps 4 and 5
    PIME_list = unlist(c(PIME_Genus_Control, PIME_Genus_T1D))
    Deseq_list = unlist(c(Deseq_Genus_Control, Deseq_Genus_T1D))
    PS = tax_glom(PS, "Genus") #[ 184 taxa and 1467 samples ]
    
    P_Control_list = PIME_Genus_Control
    P_T1D_list = PIME_Genus_T1D
    D_Control_list = Deseq_Genus_Control
    D_T1D_list = Deseq_Genus_T1D
    
    taxa_id_1 = taxa_identifiers[,c("ASV", "Genus")]
  } else {
    #Create lists of Taxa of interests from steps 4 and 5
    PIME_list = unlist(c(PIME_ASV_Control, PIME_ASV_T1D))
    Deseq_list = unlist(c(Deseq_ASV_T1D, Deseq_ASV_Control))
    taxa_id_1 = taxa_identifiers[,c("ASV", "Identifier")]
    
    
    P_Control_list = PIME_ASV_Control
    P_T1D_list = PIME_ASV_T1D
    D_Control_list = Deseq_ASV_Control
    D_T1D_list = Deseq_ASV_T1D
  }
  PS
  #####################
  #Relative Abundance and Abundance 
  ps.RA = transform_sample_counts(PS, function(x) x / sum(x) )
  
  #################################################################
  #qPCR; 16s rRNA reads/g
  qpcr_data = map[,c("ID", "copies_16s_per_gram_stool")]
  row.names(qpcr_data) = qpcr_data$ID
  qpcr_data$ID = NULL
  
  #Merge the Relative Abundance table and the qpcr table 
  otu_RA = data.frame(otu_table(ps.RA))
  otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
  row.names(otu_RA_1) = otu_RA_1$Row.names
  otu_RA_1$Row.names = NULL
  
  #Multiply all the Relabundance columns by the reads/g column
  for(i in 2:length(colnames(otu_RA_1))) {
    otu_RA_1[is.na(otu_RA_1[,i]),i] = 0
    otu_RA_1[,i] <- suppressWarnings(as.integer(round(otu_RA_1[,1] * otu_RA_1[, i])))
  }
  ps.RA.qpcr_1 = ps.RA
  otu_table(ps.RA.qpcr_1) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
  current_group = "PIME"
  current_method = "Reads"
  current_list = "HLA"
  for (current_group in c("PIME", "DESEQ")) {
    final_data = data.frame()
    if (current_group == "PIME") {
      Taxa_list = PIME_list
      Control_list = P_Control_list
      T1D_list = P_T1D_list
    } else {
      Taxa_list = Deseq_list
      Control_list = D_Control_list
      T1D_list = D_T1D_list
    }
    
    for (current_method in c("Reads", "RelativeAbundace")) {
      #Read in appropriate Phyloseq object 
      if (current_method == "Reads") {
        current_ps = ps.RA.qpcr_1 
      } else {
        current_ps = ps.RA 
      }
      #####################
      #Create data.frame
      current_Taxa = data.frame(tax_table(current_ps))
      current_Taxa$ASV = rownames(current_Taxa)
      if (current_glom != "Genus") {
        current_Taxa = merge(current_Taxa, taxa_id_1, by = "ASV")
      } 
      #Replace symbol with "." to save space 
      current_Taxa[,current_glom] = gsub( "-","\\.", current_Taxa[,current_glom])
      current_Taxa = subset(current_Taxa, current_Taxa[,current_glom] %in% Taxa_list)
      rownames(current_Taxa) = current_Taxa$ASV
      
      tax_table(current_ps) = tax_table(as.matrix(current_Taxa))
      current_Taxa$ASV = rownames(current_Taxa)
      current_Taxa = current_Taxa[,c("ASV", current_glom)]
      
      #Create table of OTU counts and map data 
      count_table = data.frame(t(otu_table(current_ps)))
      count_table$ASV = row.names(count_table)
      count_table = merge(count_table, current_Taxa)
      row.names(count_table) = count_table[,current_glom]
      count_table$ASV = NULL
      count_table[,current_glom]= NULL
      count_table = data.frame(t(count_table))
      count_table$ID = row.names(count_table)
      count_table$ID = gsub("X", "", count_table$ID)
      count_table_2 = merge(count_table, map)
      
      if (current_method == "Reads") {
        Read_Count_table = count_table_2
      } else {
        RelAbun_Count_table = count_table_2
      }
      
      Output_data = data.frame()
      Sign_data = data.frame()
      #Test significance of all metadata on Genera 
      
      for (current_variable in environmental_list) {
        pvalue_data = data.frame()
        current_data = count_table_2 
        current_data = subset(current_data, is.na(current_data[,current_variable]) == F)
        
        current_genus= Taxa_list[1]
        #If there is more than one sub-category 
        if (length(unique(current_data[,current_variable]))>1) {
          for (current_genus in Taxa_list) {
            #Run kruskal.test, extract and save p-value
            ANOVA = kruskal.test(as.formula(paste(current_genus, "~", current_variable)),
                                 data = current_data)
            Pvalue = ANOVA$p.value
            P_row = data.frame("Genus" = current_genus, "Variable" = current_variable,  "Pvalue" =  Pvalue)
            pvalue_data = rbind(pvalue_data, P_row)
          }
        }
        #Adjust pvalues and subset 
        pvalue_data$FDR = p.adjust(pvalue_data$Pvalue, "fdr")
        pvalue_data = subset(pvalue_data, pvalue_data$Pvalue <= 0.05) 
        if (nrow(pvalue_data) > 0) {
          Sign_data = rbind(Sign_data, pvalue_data)
        }
      }
      Sign_data$method = current_method
      final_data = rbind(final_data, Sign_data)
    }
    
    #Create column of Genus and significant metadata variabke 
    final_data$Combined = paste(final_data$Genus, final_data$Variable, sep = "_")
    final_table = data.frame(table(final_data$Combined))
    #Keep those present in both relative abundance and total abundance 
    final_table = subset(final_table, final_table$Freq > 1)
    final_data = subset(final_data, final_data$Combined %in% final_table$Var1)
    
    Sign_data = final_data
    #####################################
    #Graphing  
    graphing_data = data.frame()
    graphing_data_1 = data.frame()
    
    #Create long data for all significant variables 
    for (current_variable in unique(Sign_data$Variable)) {
      current_subset = subset(Sign_data,Sign_data$Variable == current_variable)
      for (current_genus in unique(current_subset$Genus)) {
        for (current_method in c("Reads", "RelAbun")) {
          if (current_method == "Reads") {
            current_count = Read_Count_table
          } else {
            current_count = RelAbun_Count_table
          }
          current_subset_1 = current_count[,c(current_variable, current_genus)]
          current_subset_1$Factor = current_variable
          current_subset_1$Genus = current_genus
          current_subset_1$Group = current_method
          colnames(current_subset_1) = c("Subfactor", "Value", "Factor", "Genus", "Group")
          current_subset_1= subset(current_subset_1, is.na(current_subset_1$Subfactor) == F)
          graphing_data = rbind(graphing_data, current_subset_1)
          current_subfactor = unique(current_subset_1$Subfactor)[1]
          for (current_subfactor in unique(current_subset_1$Subfactor)) {
            current_subset_3 = subset(current_subset_1, current_subset_1$Subfactor == current_subfactor )
            current_subset_3 = subset(current_subset_3, current_subset_3$Value > 0)
            current_mean = mean(current_subset_3$Value, na.rm = T)
            current_sd = sd(current_subset_3$Value, na.rm = T)
            current_row = data.frame(current_variable, current_genus, current_subfactor, 
                                     current_mean, current_sd, current_method)
            colnames(current_row) = c("Factor", 'Genus', "Subfactor", "Mean", "SD", "Group")
            graphing_data_1 = rbind(graphing_data_1, current_row)
          }
        }
      }
    }
    
    #Polish up column values 
    graphing_data$Factor = gsub("_", " ", graphing_data$Factor)
    graphing_data$Factor = ifelse(graphing_data$Factor == "Apartment", "Residence Type", graphing_data$Factor)
    graphing_data$Factor = ifelse(graphing_data$Factor == "Alcohol Smoking Medications Pregnancy", "Risk Events Pregnancy", graphing_data$Factor)
    graphing_data$Factor = gsub(" Pregnancy", "\nPregnancy", graphing_data$Factor)
    graphing_data$Factor = gsub(" Binned", "", graphing_data$Factor)
    graphing_data$Factor = gsub(" Breastfeeding", "\nBreastfeeding", graphing_data$Factor)
    graphing_data$Factor = gsub(" ChronicIllness", "\nChronicIllness", graphing_data$Factor)
    
    graphing_data$Factor_Total = paste(graphing_data$Factor, graphing_data$Subfactor, sep = "\n")
    graphing_data$Factor_Total = gsub("TRUE", "Present", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("FALSE", "Absent", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("3.5_times", "3-5 times", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("1.2_times", "1-2 times", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("1_3", "1-3", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("4_7", "4-7", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("8_9", "8-9", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("1_6", "1-6", graphing_data$Factor_Total)
    graphing_data$Factor_Total = gsub("7_9", "7-9", graphing_data$Factor_Total)
    
    graphing_data$Subfactor = gsub(" ", "_", graphing_data$Subfactor)
    
    graphing_data$Genus = gsub("_group", " group", graphing_data$Genus)
    graphing_data$Genus = gsub("_", "\n", graphing_data$Genus)
    
    #Manually select colors for each subcatergory 
    break_list = c("Male","Female",
                   "No","Yes",
                   "1.2_times","3.5_times","Seldom", "Never",
                   "1.2_per_week","3.5_per_week","Daily",
                   "1_3","4_7","8_9",
                   "FALSE","TRUE",
                   "Other","Flat",
                   "0","1","2","3")
    value_list = c("cornflowerblue", "darksalmon", 
                   "tomato3", "springgreen4", 
                   "gray65", "gray55", "gray42", "gray22",
                   "plum2", "mediumorchid", "mediumorchid4", 
                   "deepskyblue", "dodgerblue", "dodgerblue4",
                   "red", "limegreen", 
                   "slateblue3", "turquoise3", 
                   "coral1", "brown2", "firebrick2", "firebrick4")
    
    Control_graph = ""
    T1D_graph = ""
    for (control_T1D in c("T1D", "Control")) {
      if (control_T1D == "T1D") {
        subset_list = T1D_list
      } else {
        subset_list = Control_list
      }
      
      
      if (length(subset_list) > 0 ) {
        graphing_data_subset = subset(graphing_data, graphing_data$Genus %in% subset_list)
        if (nrow(graphing_data_subset)> 0) {
          #Create list of plots   
          plot_list = list()
          i = 1 #Create starting iteration value for plot list 
          for (current_genus in sort(unique(graphing_data_subset$Genus))) {
            current_subset = subset(graphing_data_subset, graphing_data_subset$Genus == current_genus)
            
            current_subset_RelAbun = subset(current_subset, current_subset$Group != "Reads")
            current_RelAbun_plot = ggplot(current_subset_RelAbun, aes(x = Factor_Total, y = Value, color = Subfactor)) +
              geom_boxplot() + 
              coord_cartesian(ylim = c(0, mean(current_subset_RelAbun$Value, na.rm = 2)*2)) +
              theme(legend.position="none",
                    axis.title.x = element_blank(),
                    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7)) +
              scale_color_manual(breaks = break_list, values= value_list) +
              ggtitle(current_genus) + ylab("RA")
            current_plot = current_RelAbun_plot
            plot_list[[i]] = current_plot
            i = i + 1 
          }
          if (current_group == "PIME") {
            final_plot = ggarrange(plotlist = plot_list, nrow = 3, ncol = 2); final_plot
          } else {
            final_plot = ggarrange(plotlist = plot_list); final_plot
          }
          
          final_plot = annotate_figure(final_plot, paste("More in ", control_T1D, sep = ""))
          assign(paste(control_T1D, "_graph", sep = ""), final_plot)
        }
      } 
    }
    if (typeof(T1D_graph) == "list" & typeof(Control_graph) == "list") {
      final_plot = ggarrange(T1D_graph, Control_graph, ncol = 2)
    } else if (typeof(T1D_graph) == "character" & typeof(Control_graph) == "list" ) {
      final_plot = Control_graph
    } else if (typeof(T1D_graph) == "list" & typeof(Control_graph) == "character") {
      final_plot = T1D_graph
    }
    final_plot = annotate_figure(final_plot, paste(current_group, "\nRelative Abundance"))
    #Assign and save dataframe and plot
    assign(paste(current_group, "_plot", sep = ""), final_plot)
    print(paste(current_group, "_plot", sep = ""))
    assign(paste(current_group, "_data", sep = ""), Sign_data)
    
  }
  
  #Save dataframes of: Genus, Variable, Pvaluem FDR
  write.csv(DESEQ_data, paste("./CSV_files/DESEQ_confounder_", current_glom, ".csv", sep = ""), row.names = F)
  write.csv(PIME_data, paste("./CSV_files/PIME_confounder_", current_glom, ".csv", sep = ""), row.names = F)
  
  #Save plots as Jpegs
  jpeg(paste("./Images/PIME_Confounders_", current_glom, ".jpeg", sep = ""), res = 400, height = 6000, width = 5000)
  print(PIME_plot)
  dev.off()
  
  jpeg(paste("./Images/DESEQ_Confounders_", current_glom, ".jpeg", sep = ""), res = 400, height = 4000, width = 5000)
  print(DESEQ_plot)
  dev.off()
  
  assign(paste(current_glom, "_Deseq", sep = ""), DESEQ_plot)
  assign(paste(current_glom, "_PIME", sep = ""), PIME_plot)
  
}


jpeg("./Images/Genus_Confounders.jpeg",
     res = 400, height = 4000, width = 6000)
ggarrange(Genus_PIME, Genus_Deseq, ncol = 2, 
          widths = c(2, 1))
dev.off()



