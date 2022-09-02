#####################
#Overview 
#####################
#Create Alpha Diversity Plots of Total Population and Matched Controls

#####################
#Libraries and Functions 
#####################
'%notin%' = Negate('%in%') #https://stackoverflow.com/questions/38351820/negation-of-in-in-r
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
library(ggplot2) #https://ggplot2.tidyverse.org/
library(ggpubr) #https://cran.r-project.org/web/packages/ggpubr/index.html

#####################
#Load Files 
#####################
setwd("~/Desktop/T1D_Demultiplexed/")

#From Subsetting Step 
ps.filtered = readRDS("Filters_PS.RDS") #Also works for Relative Abundance
Total_Map = read.csv("Matched_Samples.csv")
rownames(Total_Map) = Total_Map$ID
sample_data(ps.filtered) = sample_data(Total_Map)

Total_Map = data.frame(sample_data(ps.filtered)) #1223
Total_Map$Autoimmune_2_groups = ifelse(Total_Map$Autoimmune_2_groups == "Control", "All Control", "Type 1\nDiabetes")

table(Total_Map$Autoimmune_2_groups)
current_glom = "Genus"
for (current_glom in c("Genus", "ASV")) {
  if (current_glom == "Genus") {
    ps.filtered_1 = tax_glom(ps.filtered, "Genus")
  } else {
    ps.filtered_1 = ps.filtered
  }
  ##########
  #Make Reads/g phyloseq object 
  ps.RA = transform_sample_counts(ps.filtered_1, function(x) x / sum(x) )
  
  #################################################################
  #qPCR; 16s rRNA reads/g
  qpcr_data = Total_Map[,c("ID", "copies_16s_per_gram_stool")]
  row.names(qpcr_data) = qpcr_data$ID
  qpcr_data$ID = NULL
  
  #Merge the Relative Abundance table and the qpcr table 
  otu_RA = data.frame(otu_table(ps.RA))
  otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
  row.names(otu_RA_1) = otu_RA_1$Row.names
  otu_RA_1$Row.names = NULL
  
  #Multiply all the Relabundance columns by the reads/g column
  for(i in 2:length(colnames(otu_RA_1))) {
    otu_RA_1[,i] <- as.integer(otu_RA_1[,1] * otu_RA_1[, i])
  }
  
  ps.RA.qpcr = ps.RA
  otu_table(ps.RA.qpcr) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
  
  #####################
  #Total Population 
  #####################
  #Create new Phyloseq Object; using All 1519 Subjects 
  ps.current = ps.RA.qpcr
  #Rename the Cases vs Controls
  sample_data(ps.current) = sample_data(Total_Map)
  
  #Create alpha diversity boxplots 
  #https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness
  #https://www.rdocumentation.org/packages/ggpubr/versions/0.4.0/topics/stat_compare_means
  Matched_All = plot_richness(ps.current, x="Autoimmune_2_groups", measures=c("Shannon", "Observed")) +
    geom_boxplot() + stat_compare_means() + 
    ggtitle(current_glom)+
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
          axis.title.x = element_blank())
  
  assign(paste(current_glom, "plot", sep = "_"), Matched_All)
  print(paste(current_glom, "plot", sep = "_"))
}

#Save combined figure as jpeg in working directory 
#https://astrostatistics.psu.edu/su07/R/html/grDevices/html/png.html
jpeg("./Images/Alpha.jpeg", height = 3000, width = 4000, res = 400) 
ggarrange(Genus_plot, ASV_plot, nrow = 2)
dev.off()
