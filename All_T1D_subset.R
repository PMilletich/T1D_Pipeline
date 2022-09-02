setwd("~/Desktop/T1D_Demultiplexed/")
Sample_data = read.csv("T1D_Samples.csv")


subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")
Subset_1 = c("Region","Siblings_at_birth", "Apartment", "Total_Breastfeeding_Binned")
Matched_controls = subset_cases

count = c()
current_id = unique(subset_cases$ID)[1]
for (current_id in unique(subset_cases$ID)){
  current_row = subset(subset_cases, subset_cases$ID == current_id)
  current_row_subset = current_row[,Subset_1]
  current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
  
  Matched_subset = merge(subset_controls, current_row_subset)
  count = c(count, nrow(Matched_subset))
  Matched_controls = rbind(Matched_controls, Matched_subset)
}
Matched_controls = unique(Matched_controls)
table(Matched_controls$Autoimmune_2_groups)

mean(count)
sd(count)

write.csv(Matched_controls, "Matched_Samples.csv")

nrow(subset(Sample_data, 
            is.na(Sample_data$Exclusive_Breastfeeding_Binned) == F))/nrow(Sample_data)
