#####################
#Overview 
#####################
#Find Factors impacting Autoimmunity - Chisq.test
#find Factors impacting Binomial Beta Diversity - Adonis

#####################
#Libraries and Functions 
#####################
library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
'%notin%' = Negate('%in%')#https://stackoverflow.com/questions/38351820/negation-of-in-in-r

#####################
#Load Files 
#####################
setwd("~/Desktop/T1D_Demultiplexed/")
#Sample Data and Phyloseq from Previous Subsetting Step 
OG_PS = readRDS("Filters_PS.RDS")
map = data.frame(sample_data(OG_PS))

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

#################################################################
#qPCR 
#################################################################
qpcr_data = map[,c("ID", "copies_16s_per_gram_stool")]
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

ps.qpcr = ps.RA
otu_table(ps.qpcr) = otu_table(ps.qpcr, taxa_are_rows = F)

table(map$Autoimmune_2_groups)
################################
#Create list of columns to compare 
################################

HLA = c("DR13.DQ603","DR13.DQ604","DR15.DQ602","DR3.DQ2.5","DR4.DQ7",
        "DR4.DQ8","DR5.DQ7","DR7.DQ2","DR7.DQ9","DR8.DQ4","DR9.DQ9","DR1.DQ5")

Immunlogical = c("Antibiotics_during_Pregnancy","Infection_during_Pregnancy","Infection_as_newborn","Infection_w_antibiotics_1.12m",
                 "Cold.UpperRespTractInf","Otitis","Pneumonia","gastroenteritis","OtherInfection","OtherDisease",
                 "Antibiotics_Pregnancy","Corticosteroids_Pregnancy","HighBloodPressure_Meds_Pregnancy","Psychiatric_Meds_Pregnancy",
                 "PainKiller_Meds_Pregnancy","HormonePreparates_Pregnancy","Other_Meds_Pregnancy","Total_Meds_Pregnancy")

Dietary = c("Intro_CowsMilk_Binned","Intro_Gluten_Binned","Egg_first_year","Beef_first_year","Pork_first_year",
            "Exclusive_Breastfeeding_Binned","Total_Breastfeeding_Binned","Intro_formula_Binned")

Other = c("Autoimmune_2_groups","Region","Sex","Siblings_at_birth","Delivery","Apartment","Father_Only_Elementary","Mother_Elementary",
          "Father_Unemployed","Mother_Unemployed","Both_Parents_Abroad","Single_Mother","Stressful_Life_Event",
          "Mother_NoSupport_Pregnancy", "Mother_NotSafe_Pregnancy","Worry_ChronicIllness_Child",
          "Smoking_Pregnancy","Alcohol_Pregnancy_PerMonth","Risky_Alcohol_Pregnancy",
          "Alcohol_Smoking_Medications_Pregnancy","Mother_Over_35","Father_Over_40" )

################################
#### Factors impacting Autoimmune Development 
################################
#create new dataframes for Cases or Controls
AI = subset(map, map$Autoimmune_2_groups != "Control")
Control = subset(map, map$Autoimmune_2_groups == "Control")

Final_df = data.frame()
Distribution = data.frame()

for (current_list_var in c("HLA", "Diet", "Immune", "Other")) {
  #All Factors and Pvalues 
  Case_Confounders = data.frame()
  #Significant Pvalues and Distributions
  if (current_list_var == "HLA") {
    environmental_list = HLA
  } else if (current_list_var == "Diet") {
    environmental_list = Dietary
  } else if (current_list_var == "Immune") {
    environmental_list = Immunlogical
  } else if (current_list_var == "Other") {
    environmental_list = Other
  } 
  for (current_variable in environmental_list) {
    #Create a new dataframe of subjects with information on the current variable 
    map_current = subset(map, is.na(map[,current_variable]) == F)
    
    #Extract Pvalue from chisq.test
    #https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/chisq.test
    pvalue = chisq.test( map_current[,current_variable], map_current$Autoimmune_2_groups, 
                         simulate.p.value = TRUE)$p.value
    
    #Bind the pvalue and variable to master dataframe 
    current_row = data.frame("Variable" = current_variable, 
                             "Pvalue" = pvalue)
    Case_Confounders = rbind(Case_Confounders, current_row)
    
    #If the variable is significant 
    if (round(pvalue,3) <= 0.05) {
      #Create dataframes of the counts
      AI_Count = data.frame(table(AI[,current_variable])); colnames(AI_Count) = c("Var", "Count")
      NA_Count = nrow(AI)-sum(AI_Count$Count)
      NA_Count = data.frame("Var" = "NaN", "Count" = NA_Count)
      AI_Count = rbind(AI_Count, NA_Count)
      #Create dataframes of the percents 
      AI_Count$Percent = round(AI_Count$Count/nrow(AI), 3)*100
      
      #Paste values into one column 
      AI_Count$AI_Overall = paste(AI_Count$Count, " (", AI_Count$Percent, "%)", sep = "")
      #Only keep variable and pasted column 
      AI_Total = AI_Count[,c("Var", "AI_Overall")]
      
      #Repeat for controls
      Control_Count = data.frame(table(Control[,current_variable])); colnames(Control_Count) = c("Var", "Count")
      NA_Count = nrow(Control)-sum(Control_Count$Count)
      NA_Count = data.frame("Var" = "NaN", "Count" = NA_Count)
      Control_Count = rbind(Control_Count, NA_Count) 
      Control_Count$Percent = round(Control_Count$Count/nrow(Control), 3)*100
      
      Control_Count$Control_Overall = paste(Control_Count$Count, " (", Control_Count$Percent, "%)", sep = "")
      Control_Total = Control_Count[,c("Var", "Control_Overall")]
      
      #Merge the Case and Control tables 
      AI_Control = merge(Control_Total, AI_Total, all = T)
      #Create a dataframe to bind to master dataframe
      AI_Control$Factor = current_variable
      
      Distribution = rbind(Distribution, AI_Control)
    }
  }
  Case_Confounders$FDR = p.adjust(Case_Confounders$Pvalue, "fdr")
  Final_df = rbind(Final_df, Case_Confounders)
}

#Case_Confounders = Case_Confounders[order(Case_Confounders$Pvalue),]
write.csv(Final_df, "./CSV_files/Chisq_Pvalues.csv", row.names = F)
write.csv(Distribution, "./CSV_files/Chisq_Distribution.csv", row.names = F)

################################
#### Factors impacting Binomial Beta Diversity
#### Takes a long time to run, depending on number of columns comparing 
################################
current_ps = ps.qpcr #Set Phyloseq object to Relative Abun 
Binomial_df = data.frame()

for (current_list_var in c("HLA", "Diet", "Immune", "Other")) {
  print(current_list_var)
  #All Factors and Pvalues 
  Case_Confounders = data.frame()
  #Significant Pvalues and Distributions
  Distribution = data.frame()
  if (current_list_var == "HLA") {
    environmental_list = HLA
  } else if (current_list_var == "Diet") {
    environmental_list = Dietary
  } else if (current_list_var == "Immune") {
    environmental_list = Immunlogical
  } else if (current_list_var == "Other") {
    environmental_list = Other
  } 
  
  for (current_variable in environmental_list) {
    ps.RA_Variable = current_ps
    #Creaet a new sample dataframe and remove subjects lacking the current_variable
    metadata_Variable = as(sample_data(ps.RA_Variable), "data.frame")
    metadata_Variable = subset(metadata_Variable, is.na(metadata_Variable[,current_variable]) == F)
    #Load in subsetted data to PS
    sample_data(ps.RA_Variable) = sample_data(metadata_Variable)
    
    #Create a distance matrix from the phyloseq object 
    #https://joey711.github.io/phyloseq/distance.html
    dist.uf <- phyloseq::distance(ps.RA_Variable, method = "binomial",na.rm = TRUE)
    
    #Run anova between the distance matrix and the current_variable
    #https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis
    ANOVA = adonis2(as.formula(paste("dist.uf ~ ", current_variable, sep = "")),
                    data = metadata_Variable)
    
    #Extract the Pvalue
    Pvalue = ANOVA$`Pr(>F)`[1]
    Fvalue = round(ANOVA$`F`[1],3)
    R2 = round(ANOVA$R2[1],3)
    DF = ANOVA$Df[1]
    SS = round(ANOVA$SumOfSqs[1],1)
    
    #Append the Pvalue and the variable to the master dataframe 
    current_row = data.frame("Variable" = current_variable,
                             "DF" = DF, 
                             "SumofSquares" = SS, 
                             "R2" = R2, 
                             "F" = Fvalue,
                             "Pvalue" = Pvalue)
    print(current_row)
    Case_Confounders = rbind(Case_Confounders, current_row)
  }
  Case_Confounders$padj = p.adjust(Case_Confounders$Pvalue, "fdr")
  Binomial_df = rbind(Binomial_df, Case_Confounders)
}
#Save dataframe to csv 
write.csv(Binomial_df, "./CSV_files/Beta_Diversity.csv", row.names = F)
