library(phyloseq)
library(ggplot2)
library(ggpubr) #stat_compare_means
library(gridExtra) #grid.arrange
library(vegan) #Adonis
library(pairwiseAdonis)
library(stringr) #str_replace
library(dplyr) # summary_all, sample_n
library(nortest)
library(tidyr)
library(plyr)

'%notin%' = Negate('%in%')
setwd("~/Desktop/Reviewer_Edits_v2/")


###########3
#Total_breastfeeding.mo. 
Sample_data = read.csv("Celiac_sample_data_5.csv")
Matched_data = read.csv("Matched_Celiac.csv")

environmental_list = c("Celiac_DQ")

final_data = data.frame()

current_environ = "Sex"
for (current_group in c("Celiac", "Total", "Matched")) {
  if (current_group == "Matched") {
    current_data = subset(Matched_data, Matched_data$Autoimmune_Disorder == "Control")
  } else if (current_group == "Celiac") {
    current_data = subset(Sample_data, Sample_data$Autoimmune_Disorder != "Control")
  } else if (current_group == "Total") {
    current_data = subset(Sample_data, Sample_data$Autoimmune_Disorder == "Control")
  } else {
    print("Uh oh")
  }
  
  for (current_environ in environmental_list) {
    current_subset = data.frame(table(current_data[,current_environ]))
    NA_row = data.frame("Var1" = "NA", "Freq" = nrow(current_data) - sum(current_subset$Freq))
    current_subset = rbind(current_subset, NA_row)
    current_subset$Percent = round(current_subset$Freq/nrow(current_data),2)*100
    current_subset$Print = paste(current_subset$Freq, " (", current_subset$Percent, "%)", sep = "")
    current_subset$Group = current_group
    current_subset$Factor = current_environ
    current_subset = current_subset[,c("Group", "Factor", "Var1", "Print")]
    final_data= rbind(final_data, current_subset)
  }
}
