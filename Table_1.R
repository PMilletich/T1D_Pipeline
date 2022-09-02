###### 
#File import 
setwd("~/Desktop/T1D_Demultiplexed/")
Sample_data = read.csv("Matched_Samples.csv")

T1D = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
Controls_All = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

#Matching list 
Subset_1 = c("Region","Siblings_at_birth", "Apartment", "Total_Breastfeeding_Binned")

Control_Num = 2
seed_number = 10

#####################
#Sample Matching with n Controls
subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

final_subset = data.frame()

Matched = list()
Count_1 = 0
current_id = "2061"#unique(subset_df_1$ID)[1]
for (current_id in unique(subset_cases$ID)){
  current_row = subset(subset_cases, subset_cases$ID == current_id)
  current_row_subset = current_row[,Subset_1]
  current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
  Matched_subset = merge(subset_controls, current_row_subset)
  Matched = c(Matched, nrow(Matched_subset))
  # print(paste(current_id, nrow(Matched_subset)))
  
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
table(final_subset$T1D.Risk.group)

Matched_Controls = final_subset

factor_list = c("T1D.Risk.group", "Sex", 
                "Region","Siblings_at_birth","Apartment",
                "Total_Breastfeeding_Binned", 
                "DR4.DQ8", "DR3.DQ2.5", "DR15.DQ602", "DR1.DQ5")

#factor_list = c("DR4.DQ8")
current_group = "Matched"
current_factor = "T1D.Risk.group"
for (current_group in c("T1D", "Controls", "Matched")) {
  print(current_group)
  if (current_group == "T1D") {
    current_data = T1D
  } else if (current_group == "Controls") {
    current_data = Controls_All
  } else {
    current_data = Matched_Controls
  }
  
  for (current_factor in factor_list) {
    current_table = data.frame(table(current_data[,current_factor]))
    if (nrow(current_table) > 0 ) {
      current_table$Percent = round(current_table$Freq/nrow(current_data),3)*100
      current_table$Values = paste(current_table$Freq, "(", 
                                   current_table$Percent, "%)", sep = "")
      print(current_factor)
      print(current_table[,c("Var1", "Values")])
    }
  }
}

summary(T1D$Age.y._at_diagnosis)
