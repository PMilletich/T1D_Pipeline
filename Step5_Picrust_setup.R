#Libraries 
library(seqinr)

PS = readRDS("Filters_PS.RDS")
Sample_data = read.csv("Matched_Samples.csv")
taxa_id = read.csv("Taxa_identifiers_v2.csv")
taxa_id = subset(taxa_id, taxa_id$Genus != "NA")
taxa_id = taxa_id[,c("ASV", "Identifier")]
rownames(Sample_data) = Sample_data$ID
T1D = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")

PS_T1D = PS
sample_data(PS_T1D) = sample_data(Sample_data)

###################
#Normalize 
#Total vs Relative Abundance 
###################
ps.RA = transform_sample_counts(PS_T1D, function(x) x / sum(x) )

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

############
#count 
############
T1D_OTU = data.frame(otu_table(ps.RA.qpcr_1))
#1965 --> 746
T1D_OTU = T1D_OTU[, colSums(T1D_OTU != 0) > 0]
T1D_OTU = data.frame(t(T1D_OTU))
T1D_OTU$ASV = rownames(T1D_OTU)
T1D_OTU = merge(T1D_OTU, taxa_id)
rownames(T1D_OTU) = T1D_OTU$Identifier
#T1D_OTU = subset(T1D_OTU, grepl("NA_", T1D_OTU$Identifier) == F)
T1D_OTU$ASV = NULL
T1D_OTU$Identifier = NULL
rownames(T1D_OTU) = gsub("X", "", rownames(T1D_OTU))

write.table(T1D_OTU, "T1D_count_TotalAbun.tsv", sep = "\t")

############
#Taxa
############
T1D_Tax = data.frame(tax_table(ps.RA.qpcr_1))
T1D_Tax$ASV = rownames(T1D_Tax)
T1D_Tax = merge(T1D_Tax, taxa_id)
rownames(T1D_Tax) = T1D_Tax$Identifier
T1D_Tax = T1D_Tax[rownames(T1D_OTU),]
T1D_Tax = T1D_Tax[,c("Identifier", "ASV")] 
#T1D_Tax = subset(T1D_Tax, grepl("NA_", T1D_Tax$Identifier) == F)
colnames(T1D_Tax) = c("name", "seq")
T1D_Tax$name = paste(">", T1D_Tax$name, sep = "")
T1D_Tax = subset(T1D_Tax, T1D_Tax$name != ">NA")
#Save as FASTA 
D <- do.call(rbind, lapply(seq(nrow(T1D_Tax)), function(i) t(T1D_Tax[i, ])))
write.table(D, file = "T1D_TotalAbun.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

