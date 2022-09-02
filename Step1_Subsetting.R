library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

'%notin%' = Negate('%in%')

#Read in original Phyloseq data and original metadata 
OG_PS = readRDS("phyloseq.13.Final.RDS") 
Start = read.csv("Regions_MetaData_2022_Neuro.csv")
#Combine phyloseq object and sample data 
rownames(Start) = Start$ID  
sample_data(OG_PS) = sample_data(Start)

###########################################################################
#Filter out those with less than 1000 raw reads 
##########################################################################
#Create an OTU table of raw reads, sum counts, filter, and reformat phyloseq object with new OTU
otu_21 = data.frame(otu_table(OG_PS))
otu_21$Counts = rowSums(otu_21)
summary(otu_21$Counts)
otu_212 = subset(otu_21, otu_21$Counts > 1000)
otu_table(OG_PS) = otu_table(otu_212, taxa_are_rows = F)

#Remove samples without qPCR data, and create new phyloseq object 
With_Gut_only_1000 = data.frame(sample_data(OG_PS))
With_Gut_only_1000 = subset(With_Gut_only_1000, is.na(With_Gut_only_1000$copies_16s_per_gram_stool) == F)
sample_data(OG_PS) = sample_data(With_Gut_only_1000)

###########################################################################
 #Remove Autoimmune 
table(With_Gut_only_1000$Autoimmune_Disorder)
###########################################################################
#Save those that are either T1D, or controls 
With_Gut_T1D_only = subset(With_Gut_only_1000, grepl("T1D", With_Gut_only_1000$Autoimmune_Disorder) == T |
                                With_Gut_only_1000$Autoimmune_Disorder == "Control")
#Save Sample data to phyloseq object
rownames(With_Gut_T1D_only) = With_Gut_T1D_only$ID
sample_data(OG_PS) = sample_data(With_Gut_T1D_only)
#Check to make sure there are no errors 
table(With_Gut_T1D_only$Autoimmune_Disorder)

###########################################################################
#Remove Neurological  
###########################################################################
#Save those that are either T1D, or controls 
Neurological = subset(With_Gut_T1D_only, With_Gut_T1D_only$Neurological == "Control" | 
                        With_Gut_T1D_only$Autoimmune_Disorder == "T1D")
#Save Sample data to phyloseq object
rownames(Neurological) = Neurological$ID
sample_data(OG_PS) = sample_data(Neurological)

###########################################################################
#Filter and Save 
###########################################################################
#Filter low prevalent ASVS 
ps.filtered = filter_taxa(OG_PS, function(x) sum(x > 5) >= 5, TRUE)

#Save Final metadate CSV and phyloseq object in working directory 
write.csv(Neurological, "T1D_Samples.csv", row.names = F)
saveRDS(ps.filtered, "Filters_PS.RDS")
