
# ***************** DATA ANALYSES **************************** ----------------------------
# Manuscript:   MICROBIOME ASSOCIATED MHAPLA POPULATIONS WITH PARASITIC VARIABILITY IN FIELD AND GREENHOUSE SOILS
# Authors:      ISAAC LARTEY, GREGORY BONITO, GIAN MARIA NICCOLO BENUCCI, TERENCE MARSH, AND HADDISH MELAKEBERHAN.
# Affiliation:  Michigan State University
# Journal:      TBD
# Date:         N/A
# ******************************************************************** ----------------------------

#****************WORKING ENVIRONMENT SETUP***************************-----------------------------

options(scipen = 999) #to use decimals
options(max.print=100000000) # to print more lines on the display
options(verbose=TRUE)

setwd("C:/Users/ISAAC/OneDrive - Michigan State University/WORK/phd/5 MANUSCRIPTS/5 Individual Mhapla microbiome")

#BiocManager::install("phyloseq")
library(phyloseq); ##packageVersion("phyloseq")
library(ape)  ##to read .tree file
library(colorspace)
library(stringi)
library(rhdf5)
library(zlibbioc)
library(S4Vectors)
library(Biostrings)
library(yaml)
library(colorspace)
library(ggplot2)
library(indicspecies)
library(vegan)
library(rlang) #phyloseq dependency
library(devtools)

#####  IMPORTING PROKARYOTES DATA INTO PHYLOSEQ AFTER DENOISING IN QIIME2 (DEBLUR)  #######
library(qiime2R)#Import .qza qiime2 files into R
library(Biostrings)#Allow the reading of fasta file

#Combine OTU table, tree file and metadata into a phyloseq file
physeq<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza", metadata = "meta_16s.tsv")
physeq

#Prepare constax_taxonomy physeq
library(stringr)
#constax_tax <- read.delim("constax_taxonomy.txt", header=TRUE, row.names=1)
#str_remove(constax_tax,"[_1]")
#constax_tax
#constax_tax <- tax_table(as.matrix(constax_tax))

silva_tax <- read.delim("silva_taxonomy.txt", header=TRUE, row.names=1)
#str_remove(silva_tax,"[_1]")
#constax_tax
silva_tax <- tax_table(as.matrix(silva_tax))




#Prepare to add sequences to physeq
sequences <- readDNAStringSet("dna-sequences.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)


#Combine phyloseq object and other phyloseq components
physeq <- merge_phyloseq(physeq, sequences, silva_tax)


#install.packages('devtools')
#devtools::install_github("mhahsler/rBLAST")


##########   DATA OPERATIONS    ###############

#Remove un-needed samples
#otu_table(physeq) <- subset(otu_table(physeq), select = -c(baGCRep1, baGCRep2, baGCRep4, baGCRep5, baFCRep1, baFCRep2, baFCRep3, baFCRep4, baFCRep5))
            
# Checking for non-bacteria taxonomies -------there was no archaea present
unique(as.data.frame(tax_table(physeq))$Kingdom) 
unique(as.data.frame(tax_table(physeq_greenhouse))$Phylum) 
unique(as.data.frame(tax_table(physeq))$Class) 
unique(as.data.frame(tax_table(physeq))$Order)
unique(as.data.frame(tax_table(physeq))$Family)
unique(as.data.frame(tax_table(physeq_greenhouse))$Genus)
unique(as.data.frame(tax_table(physeq))$Species)

#Filtering out for non-bacteria taxonomies
#physeq <- subset_taxa(physeq, Kingdom!="")
#physeq <- subset_taxa(physeq, Kingdom!="Eukaryota_1")
physeq <- subset_taxa(physeq, Order!="Chloroplast")
physeq <- subset_taxa(physeq, Family!="Mitochondria")
tax_table(physeq) 


# relabeling taxonomies --------------------------------------------------------------------------
#tax_table(physeq)[, "Kingdom"] <- gsub("Rank_1", "", tax_table(physeq)[, "Kingdom"])
#tax_table(physeq)[, "Phylum"] <- gsub("Rank_2", "", tax_table(physeq)[, "Phylum"])
#tax_table(physeq)[, "Class"] <- gsub("Rank_3", "", tax_table(physeq)[, "Class"])
#tax_table(physeq)[, "Order"] <- gsub("Rank_4", "", tax_table(physeq)[, "Order"])
#tax_table(physeq)[, "Family"] <- gsub("Rank_5", "", tax_table(physeq)[, "Family"])
#tax_table(physeq)[, "Genus"] <- gsub("Rank_6", "", tax_table(physeq)[, "Genus"])
#tax_table(physeq)[, "Species"] <- gsub("Rank_7", "", tax_table(physeq)[, "Species"])


# relabeling Incertae_sedis taxonomies to unclassified --------------------------------------------------------------------------
#physeq <- subset_taxa(physeq, Phylum!="Auriculariales_Incertae_sedis")
#physeq <- subset_taxa(physeq, Class!="Botryosphaeriales_Incertae_sedis")
#physeq <- subset_taxa(physeq, Order!="Diaporthales_Incertae_sedis")
#physeq <- subset_taxa(physeq, Family!="Diaporthales_Incertae_sedis")
physeq <- subset_taxa(physeq, Genus!="uncultured")
#physeq <- subset_taxa(physeq, Species!="Diaporthales_Incertae_sedis")


# Name unlabeled taxa as unclassified
tax_table(physeq)[tax_table(physeq)=="#VALUE!"]<- NA
tax_table(physeq)[tax_table(physeq)==""]<- NA
tax_table(physeq)[tax_table(physeq)=="uncultured"]<- NA
tax_table(physeq)[is.na(tax_table(physeq))]<-"Unclassified"

?is.name

###remove contaminants--------------------------------------------------
# using decontnam package to remove contaminants in negative controls by comparing
# distributions to those found in normal samples
library(devtools)
library(processx)
#devtools::install_github("benjjneb/decontam", force = TRUE)
#devtools::install_github("bryandmartin/corncob")
library(decontam)

#BiocManager::install("decontam")

sample_data(physeq)
write.csv(sample_data(physeq),file ="sample_check.csv")

# check library size distribution
df_physeq <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df_physeq$LibrarySize <- sample_sums(physeq)
df_physeq <- df_physeq[order(df_physeq$LibrarySize),]
df_physeq$Index <- seq(nrow(df_physeq))
write.csv(df_physeq, file = "rank_sums.csv")
ggplot(data=df_physeq, aes(x=Index, y=LibrarySize, color=Sample)) + geom_point()

# filter by prevelance 
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample == "Control"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample == "Worm", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                           contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Control)") + ylab("Prevalence (Worm)")

# remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
ps.noncontam # with contaminants removed
otu_table(ps.noncontam)

ps.noncontam

physeq_filt <- subset_samples(ps.noncontam, Sample=="Worm") #Subsetting samples to removes negative control

physeq_field <- subset_samples(physeq_filt, Origin=="Field") #Subsetting on field populations
physeq_greenhouse <- subset_samples(physeq_filt, Origin=="Greenhouse") #Subsetting on greenhouse populations


#Renaming genera based on function from literature 
library(plyr)

#tax_table(physeq)[tax_table(physeq)=="uncultured"]<- NA


physeq_filt_func <- physeq_filt
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Acidibacter"]<-"Iron_reducing"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Acidothermus"]<-"Anti_fungi"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Unclassified"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Actinophytocola"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Actinospica"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Rhizobium"]<-"Nitrogen_fixer"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Amycolatopsis"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Bradyrhizobium"]<-"Nitrogen_fixer"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Brevundimonas"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Bryobacter"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Paraburkholderia"]<-"Anti_fungi"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Candidatus Phytoplasma"]<-"Plant_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Candidatus Udaeobacter"]<-"Antibiotic_resistant"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Catenulispora"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Caulobacter"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Cellvibrio"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Chitinophaga"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Chryseolinea"]<-"Anti_fungi"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Clostridium sensu stricto 1"]<-"Plant_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Devosia"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Dokdonella"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Dongia"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Duganella"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Dyella"]<-"Plant_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Flavobacterium"]<-"Animal_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Fluviicola"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Frankia"]<-"Nitrogen_fixer"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Gemmata"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Haliangium"]<-"Anti_fungi"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Halomonas"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Holdemanella"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Hyphomicrobium"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Inquilinus"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Kibdelosporangium"]<-"Anti_bacteria"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Kribbella"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Labrys"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Lechevalieria"]<-"Anti_bacteria"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Limnohabitans"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Massilia"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Mesorhizobium"]<-"Nitrogen_fixer"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Mucilaginibacter"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Mycobacterium"]<-"Animal_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Mycoplasma"]<-"Plant_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Niastella"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Nocardia"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Nocardioides"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Novosphingobium"]<-"Enhanced_nematode_parasitism"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Ohtaekwangia"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Paenibacillus"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Pedomicrobium"]<-"Nematicidal"
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Pedobacter"]<-"Nematicidal"
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Phenylobacterium"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Polaromonas"]<-"Soybean_cyst_associated"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Pseudomonas"]<-"Root_knot_nematode_associated"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Pseudonocardia"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Rheinheimera"]<-"Anti_bacteria"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Rhizobacter"]<-"Plant_pathogenic"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Rhodomicrobium"]<-"Nitrogen_fixer"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Rhodoplanes"]<-"Root_knot_nematode_associated"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Roseiarcus"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="SM1A02"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Solirubrobacter"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Sphingobium"]<-"Polysaccharide_degrader"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Sphingomonas"]<-"Plant_growth_promoter"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Steroidobacter"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Streptomyces"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="uncultured"]<-"Other"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Variovorax"]<-"Suppressive_soils"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Vibrio"]<-"Nematicidal"	
tax_table(physeq_filt_func)[tax_table(physeq_filt_func)=="Xanthomonas"]<-"Nematicidal"	

physeq_filt_func@tax_table

physeq_field_func <- subset_samples(physeq_filt_func, Origin=="Field") #Subsetting on field populations
physeq_greenhouse_func <- subset_samples(physeq_filt_func, Origin=="Greenhouse") #Subsetting on greenhouse populations


#function_list <- tax_table(physeq_filt)[,"Mfunction"]
new_function_list <- revalue(function_list, c( "Acidibacter"="Iron_reducing",		
           "Acidothermus"="Anti_fungi",		
           "Unclassified"="Other",		
           "Actinophytocola"="Suppressive_soils",		
           "Actinospica"="Suppressive_soils ",		
           "Rhizobium"="Nitrogen_fixer",		
           "Amycolatopsis"="Plant_growth_promoter",		
           "Bradyrhizobium"="Nitrogen_fixer",		
           "Brevundimonas"="Plant_growth_promoter",		
           "Bryobacter"="Other",		
           "Paraburkholderia"="Anti_fungi",		
           "Candidatus Phytoplasma"="Plant_pathogenic",		
           "Candidatus Udaeobacter"="Antibiotic_resistant",		
           "Catenulispora"="Other",		
           "Caulobacter"="Plant_growth_promoter",		
           "Cellvibrio"="Nematicidal",		
           "Chitinophaga"="Nematicidal",		
           "Chryseolinea"="Anti_fungi",		
           "Clostridium sensu stricto 1"="Plant_pathogenic",		
           "Devosia"="Nematicidal",		
           "Dokdonella"="Other",		
           "Dongia"="Other",		
           "Duganella"="Nematicidal",		
           "Dyella"="Plant_pathogenic",		
           "Flavobacterium"="Insect_pathogenic",		
           "Fluviicola"="Other",		
           "Frankia"="Nitrogen_fixer",		
           "Gemmata"="Nematicidal",		
           "Haliangium"="Anti_fungi",		
           "Halomonas"="Plant_growth_promoter",		
           "Holdemanella"="Other",		
           "Hyphomicrobium"="Other",		
           "Inquilinus"="Other",		
           "Kibdelosporangium"="Anti_bacteria",		
           "Kribbella"="Suppressive_soils ",		
           "Labrys"="Other",		
           "Lechevalieria"="Anti_bacteria",		
           "Limnohabitans"="Other",		
           "Massilia"="Suppressive_soils ",		
           "Mesorhizobium"="Nitrogen_fixer",		
           "Mucilaginibacter"="Plant_growth_promoter",		
           "Mycobacterium"="Animal_pathogenic",		
           "Mycoplasma"="Plant_pathogenic",		
           "Niastella"="Suppressive_soils ",		
           "Nocardia"="Plant_growth_promoter",		
           "Nocardioides"="Nematicidal",		
           "Novosphingobium"="Enhanced_nematode_parasitism",		
           "Ohtaekwangia"="Suppressive_soils",		
           "Paenibacillus"="Plant_growth_promoter",		
           "Pedomicrobium"="Nematicidal",		
           "Phenylobacterium"="Nematicidal",		
           "Polaromonas"="Soybean_cyst_associated",		
           "Pseudomonas"="Root_knot_nematode_associated",		
           "Pseudonocardia"="Other",		
           "Rheinheimera"="Anti_bacteria",		
           "Rhizobacter"="Plant_pathogenic",		
           "Rhodomicrobium"="Nitrogen_fixer",		
           "Rhodoplanes"="Root_knot_nematode_associated",		
           "Roseiarcus"="Other",		
           "SM1A02"="Other",		
           "Solirubrobacter"="Nematicidal",		
           "Sphingobium"="Polysaccharide_degrader",		
           "Sphingomonas"="Plant_growth_promoter",		
           "Steroidobacter"="Suppressive_soils",		
           "Streptomyces"="Nematicidal",		
           "Variovorax"="Suppressive_soils",		
           "Vibrio"="Nematicidal",		
           "Xanthomonas"="Nematicidal"		
))

#list_mat <- tax_table(taxonomy)
#list_mat <- tax_table(new_function_list)
#function_list_physeq <- merge_phyloseq(list_mat, sample_data(physeq_filt), otu_table(physeq_filt, phy_tree(physeq_filt))) #new physoleq object for functional analysis

#physeq_field_func <- subset_samples(function_list_physeq, Origin=="Field") #Subsetting on field populations
#physeq_greenhouse_func <- subset_samples(function_list_physeq, Origin=="Greenhouse") #Subsetting on greenhouse populations

#write.csv(new_function_list, "filtlist_mat.csv")

#taxonomy = read.csv("filtlist_mat.csv", sep=",", row.names=1)
#taxonomy = as.matrix(taxonomy)

# Filtering out OTUs < than 5 reads ------------------------------------------------------------
#biom_16S_uparse -> biom_16S_uparse_filt
#otu_table(biom_16S_uparse_filt) <- otu_table(biom_16S_uparse_filt)[which(rowSums(otu_table(biom_16S_uparse_filt)) >= 5),] ### PCR Errors 
#biom_16S_uparse_filt


###################  DATA SUMMARY   #########################
##################################################################################
#setting working directory
#wdir <- "C:/Users/ISAAC/Google Drive/WORK/phd/5 MANUSCRIPTS/3 Soil Microbiome/Analysis/2018NRKNPV"
output <-(paste( wdir,"/Output/", sep = ""))
source("Alpha_div_function.r")
source("miseqR.R")


#BiocManager::install("microbiome")
library(phyloseq)
library(ggplot2)## for plotting
library(magrittr)
library(ggpubr)##for combining the plots
library(vegan)##for community ecology based codes
library(limma)
library(edgeR)
library(Hmisc)
library(igraph)
library(labdsv)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(microbiome) #to summarize phyloseq
library(knitr)
library(fantaxtic) #used for creating top taxa relative abundance
library(data.table) ## Creating data table for relative abundance


#Burkholderia
#install.packages("remotes")
#remotes::install_github("DanielSprockett/reltools")
library(reltools) # use it to output a fasta file from phyloseq
Burkholderia <- subset_taxa(physeq_filt, Genus=="Burkholderia") #subsetting Burkholderia
save_fasta(ps = Burkholderia, file = "Burkholderia.fasta", rank = "Genus") # creating a fasta




#*************************** ANOSIM BETWEEN FIELD AND GREENHOUSE POPULATIONS ***************

library(vegan)
meta_anosim <- data.frame(sample_data(physeq_filt))
otu_anosim <- as.matrix(t(otu_table(physeq_filt)))


anosim_stat = anosim(otu_anosim, meta_anosim$Origin, distance = "bray", permutations = 9999)
anosim_stat








#******************************* FIELD MHAPLA POPULATIONS *********************************#


## RAW ABUNDANCE ###

## PHYLUM
physeq_bacRm = merge_samples(physeq_field, "Population2")
sample_data(physeq_bacRm)$ Field <- levels(sample_data(physeq_filt)$Population2)

pm_Phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate taxa at order level
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Phylum)           # Sort data frame alphabetically by Phylum
pm_Phylum

write.table(pm_Phylum, "untransformed_field_16S_phyla.txt")






## GENERA
physeq_bacRm = merge_samples(physeq_field, "Population2")
sample_data(physeq_bacRm)$ Field <- levels(sample_data(physeq_filt)$Population2)

pm_Genus <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus
pm_Genus

write.table(pm_Genus, "untransformed_16S_field_genera.txt")





                    
### RELATIVE ABUNDANCE GENERA ###

launch_evo_palette()
palette_box()

#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(physeq_field, "Population")
sample_data(physeq_bacRm)$Population <- levels(sample_data(physeq_field)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)

pm_phylum.field <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus



Abundance.field <- data.table(pm_phylum.field)
Abundance.field [(Abundance <= 0.00)]#Abundance.field [(Abundance <= 0.00), Genus:= "Other"]

#positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
positions <- c("Field13","Field8","Field2","Field4","Field5","Field6","Field10","Field14","Field15")
#positions <- c("1","2","3","4","5")

physeq_field@tax_table


p <- ggplot(data=Abundance.field, aes(x=Sample, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("Acidibacter" = "#652926",
                               "Acidothermus" ="steelblue", 
                               "Actinophytocola" ="#C84248", 
                               "Actinospica" ="green",
                               "Amycolatopsis" ="#E5D8BD", 
                               "Bryobacter" ="#673770", 
                               "Bosea" ="#D14285", 
                               "Bradyrhizobium" ="#952926", 
                               "Brevundimonas" ="Blue",
                               "Candidatus Phytoplasma" = "#673770",
                               "Candidatus Udaeobacter" = "#AD6F3B", 
                               "Catenulispora" ="#CBD588", 
                               "Caulobacter" ="#5F7FC7", 
                               "Cellvibrio" ="orange" ,
                               "Chitinophaga" ="#DA5724",
                               "Chryseolinea" = "#508578", 
                               "Clostridium sensu stricto 1" ="#CD9BCD", 
                               "Devosia" ="tan", 
                               "Dokdonella" ="magenta1",
                               "Dongia" ="gray52" ,
                               "Duganella" ="darkorange4" ,
                               "Dyella" ="lightsalmon1" ,
                               "Flavobacterium" ="palevioletred3" ,
                               "Fluviicola" ="olivedrab1" ,
                               "Frankia" ="cyan1", 
                               "Gemmata" ="darkviolet" ,
                               "Halomonas" ="lavender", 
                               "Haliangium" ="greenyellow", 
                               "Holdemanella" ="sienna1",
                               "Hyphomicrobium" ="deepskyblue1" ,
                               "Inquilinus" ="honeydew2",
                               "Kibdelosporangium" = "red", 
                               "Kribbella" ="aquamarine4" ,
                               "Labrys" ="gold" ,
                               "Lechevalieria" ="plum3" ,
                               "Limnohabitans" ="peachpuff3", 
                               "Massilia" ="turquoise3" ,
                               "Mesorhizobium" ="blueviolet" ,
                               "Mucilaginibacter" ="green4", 
                               "Mycobacterium" ="palegreen2",
                               "Mycoplasma" ="#CCCCCC", 
                               "Niastella" ="hotpink3" ,
                               "Nocardia" ="gray38" ,
                               "Nocardioides" ="black" ,
                               "Novosphingobium" ="midnightblue" ,
                               "Ohtaekwangia" ="khaki3" ,
                               "Paenibacillus" ="yellow4" ,
                               "Paraburkholderia" ="pink",
                               "Pedomicrobium" ="cadetblue",
                               "Pedobacter"="purple4",
                               "Phenylobacterium" ="rosybrown3" ,
                               "Polaromonas" ="gray69" ,
                               "Pseudomonas" ="paleturquoise4",
                               "Pseudonocardia" = "limegreen" ,
                               "Rheinheimera" ="lemonchiffon3",
                               "Rhizobacter" = "deeppink1" ,
                               "Rhizobium" ="#AD6F3B",
                               "Rhodomicrobium" ="plum4" ,
                               "Rhodoplanes" ="darkolivegreen1" ,
                               "Roseiarcus" ="darkorchid1" ,
                               "SM1A02" ="rosybrown4", 
                               "Solirubrobacter" ="cornflowerblue", 
                               "Sphingobium" ="lightgoldenrodyellow", 
                               "Sphingomonas" ="springgreen1" ,
                               "Steroidobacter" ="moccasin",
                               "Streptomyces" ="lavenderblush", 
                               "Unclassified" =	"darksalmon",
                               "Variovorax" ="navyblue" ,
                               "Vibrio" ="yellow" ,
                               "Xanthomonas" ="#D9D9D9" 
                               )) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=25)) + scale_x_discrete(limits = positions) + theme_gray()











### RELATIVE ABUNDANCE PHYLA ###

physeq_bacRm = merge_samples(physeq_field, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(physeq_field)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)

pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate taxa at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Phylum)           # Sort data frame alphabetically by Phylum


dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00)] #dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Phylum:= "Other"]
write.table(dat_Soil_Pre_Abs_bp,"Phylum_rel_abun_field.csv") #output the rel. abu table for further analysis in excel


#positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
positions <- c("Field13","Field8","Field2","Field4","Field5","Field6","Field10","Field14","Field15")

p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Phylum))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("Acidobacteriota" = "#7FC97F", 
                               "Actinobacteriota" = "#BEAED4", 
                               "Bacteroidota" = "#FDC086",
                               "Chloroflexi" = "#386CB0",
                               "Firmicutes" = "#F0027F", 
                               "Myxococcota" = "#BF5B17", 
                               "Planctomycetota" = "#666666", 
                               "Proteobacteria" = "#1B9E77", 
                               "Verrucomicrobiota" = "#FF7F00" 
                                                              )) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) + scale_x_discrete(limits = positions) + theme_gray()











#RELATIVE ABUNDANCE BY FUNCTION


#launch_evo_palette()
#palette_box()

physeq_bacRm = merge_samples(physeq_field_func, "Population")
sample_data(physeq_bacRm)
#sample_data(physeq_bacRm)$Population <- levels(sample_data(physeq_field_func)$Population)
#physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
#View(physeq_bacRm)

pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Phylum


dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00)] #dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Phylum:= "Other"]

#positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
positions <- c("Field13","Field8","Field2","Field4","Field5","Field6","Field10","Field14","Field15")


p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) + 
  scale_fill_manual(values = c("Animal_pathogenic" ="deeppink1",
                               "Anti_bacteria" = "#652926",
                               "Anti_fungi" ="steelblue",
                               "Antibiotic resistant" ="#E5D8BD",
                               "Enhanced_nematode_parasitism" ="#C84248", 
                               "Iron_reducing" ="#5F7FC7", 
                               "Nematicidal" ="Blue",
                               "Nitrogen_fixer" ="green",
                               "Other" = "#673770", 
                               "Plant_growth_promoter" ="#DA5724",
                               "Plant_pathogenic" = "#AD6F3B",
                               "Polysaccharide_degrader" ="#673770", 
                               "Root_knot_nematode_associated" ="#CBD588",
                               "Soybean_cyst_associated" ="#D14285", 
                               "Suppressive_soils" ="orange"
  )) +
  theme(legend.position="right") + scale_x_discrete(limits = positions) + theme_gray() + guides(fill=guide_legend(nrow=15))






#Splitting populations to identify the presence/absence of genera-----------------------------------------------------------------
Abundance.field.split <- split(Abundance.field, Abundance.field$Sample)

#Population 13
Abundance.field.split.13 <- Abundance.field.split$Field13
Abundance.field.split.13.genera <- unique(Abundance.field.split.13$Genus)

#Population 8
Abundance.field.split.8 <- Abundance.field.split$Field8
Abundance.field.split.8.genera <- unique(Abundance.field.split.8$Genus)

#Population 2
Abundance.field.split.2 <- Abundance.field.split$Field2
Abundance.field.split.2.genera <- unique(Abundance.field.split.2$Genus)

#Population 4
Abundance.field.split.4 <- Abundance.field.split$Field4
Abundance.field.split.4.genera <- unique(Abundance.field.split.4$Genus)

#Population 5
Abundance.field.split.5 <- Abundance.field.split$Field5
Abundance.field.split.5.genera <- unique(Abundance.field.split.5$Genus)

#Population 6
Abundance.field.split.6 <- Abundance.field.split$Field6
Abundance.field.split.6.genera <- unique(Abundance.field.split.6$Genus)

#Population 10
Abundance.field.split.10 <- Abundance.field.split$Field10
Abundance.field.split.10.genera <- unique(Abundance.field.split.10$Genus)

#Population 14
Abundance.field.split.14 <- Abundance.field.split$Field14
Abundance.field.split.14.genera <- unique(Abundance.field.split.14$Genus)

#Population 15
Abundance.field.split.15 <- Abundance.field.split$Field15
Abundance.field.split.15.genera <- unique(Abundance.field.split.15$Genus)
Abundance.field.split.15.genera <- unique(Abundance.field.split.15$Genus)

#All genera in field populations
field.genera <- unique(Abundance.field$Genus)

# Finding maximum length
max_ln1 <- max(c(length(Abundance.field.split.13.genera), length(Abundance.field.split.8.genera)))
max_ln2 <- max(c(length(Abundance.field.split.2.genera), length(Abundance.field.split.4.genera)))
max_ln3 <- max(c(length(Abundance.field.split.5.genera), length(Abundance.field.split.6.genera)))
max_ln4 <- max(c(length(Abundance.field.split.10.genera), length(Abundance.field.split.14.genera)))
max_ln5 <- max(c(length(Abundance.field.split.15.genera), length(field.genera)))

max_ln<-max(max_ln1,max_ln2,max_ln3,max_ln4,max_ln5)

# Merge all populations
Pop.Field.13 = c(Abundance.field.split.13.genera,rep(NA, max_ln - length(Abundance.field.split.13.genera)))
Pop.Field.8 = c(Abundance.field.split.8.genera,rep(NA, max_ln - length(Abundance.field.split.8.genera)))
Pop.Field.2 = c(Abundance.field.split.2.genera,rep(NA, max_ln - length(Abundance.field.split.2.genera)))
Pop.Field.4 = c(Abundance.field.split.4.genera,rep(NA, max_ln - length(Abundance.field.split.4.genera)))
Pop.Field.5 = c(Abundance.field.split.5.genera,rep(NA, max_ln - length(Abundance.field.split.5.genera)))
Pop.Field.6 = c(Abundance.field.split.6.genera,rep(NA, max_ln - length(Abundance.field.split.6.genera)))
Pop.Field.10 = c(Abundance.field.split.10.genera,rep(NA, max_ln - length(Abundance.field.split.10.genera)))
Pop.Field.14 = c(Abundance.field.split.14.genera,rep(NA, max_ln - length(Abundance.field.split.14.genera)))
Pop.Field.15 = c(Abundance.field.split.15.genera,rep(NA, max_ln - length(Abundance.field.split.15.genera)))
All.field.genera = c(field.genera,rep(NA, max_ln - length(field.genera)))


Combined.field.genera <- cbind(All.field.genera, Pop.Field.13, Pop.Field.8, Pop.Field.2, Pop.Field.4, Pop.Field.5, Pop.Field.6, Pop.Field.10, Pop.Field.14, Pop.Field.15)

write.csv(Combined.field.genera, file="All.field.genus2_silva.csv")





#Mollicutes------------------------------------------------------------------------------------------------------------------
Mollicutes.field <- subset_taxa(physeq_field, Class=="Mollicutes")  ### select greenhouse or fields
Mollicutes.greenhouse <- subset_taxa(physeq_greenhouse, Class=="Mollicutes")
Mollicutes <- subset_taxa(physeq_filt, Class=="Mollicutes")

write.table(Mollicutes@refseq, "Mollicutes.txt")


aa<-refseq(Mollicutes.field)

tax_table(Mollicutes.field)


#Relative abundance fields otu table to 100%
#install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()
palette_box()


physeq_bacRm = merge_samples(Mollicutes, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(Mollicutes)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)

pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus


dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Genus:= "Other"]

#positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
positions <- c("Field13","Field8","Field2","Field4","Field5","Field6","Field10","Field14","Field15")

p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C",
                               "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" ,"#FBB4AE", "#B3CDE3" ,"#CCEBC5", "#DECBE4", "#FED9A6",
                               "#FFFFCC" ,"#E5D8BD", "#FDDAEC", "#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
                               "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999" ,"#66C2A5",
                               "#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072",
                               "#80B1D3", "#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F")) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) + scale_x_discrete(limits = positions) + theme_gray()


#install.packages("remotes")
#remotes::install_github("mhahsler/rBLAST")

#install.packages('devtools')
#devtools::install_github("mhahsler/rBLAST")
#devtools::install_bioc("Biostrings")





###    ALPHA DIVERSITY ###

#load libraries
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(vegan)
library(ggplot2)

#plot alpha diversity
ps1 <- prune_taxa(taxa_sums(physeq_field) > 0, physeq_field) ## use unfiltered phyloseq object for alpha diversity
alpha_diversity <- alpha(ps1, index = "all")
head(alpha_diversity)


#Prepare data for visualisation
ps1.meta <- meta(ps1)
head(ps1.meta)

#Add the diversity table to metadata
ps1.meta$diversity_observed  <- alpha_diversity$observed 
ps1.meta$diversity_shannon <- alpha_diversity$diversity_shannon

#Observed diversity/Richness
ggplot(ps1.meta, aes(x = Population, y = diversity_observed)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.
kruskal.test(diversity_observed ~ Population, data=ps1.meta) #p-value = 0.2566


#Shannon diversity
ggplot(ps1.meta, aes(x = Population, y = diversity_shannon)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.
kruskal.test(diversity_shannon ~ Population, data=ps1.meta) #p-value = 0.2161






###   BETA DIVERSITY   ####

# any sample with less than 5 reads for a particular otu will be placed to 0
physeq_bac_filter <- physeq_field
#otu_table(physeq_bac_filter)[otu_table(physeq_bac_filter) <= 5] < 0- 0 ### tag switching
#otu_table(biom_16s_qc) <- otu_table(biom_16s_qc)[rowSums(otu_table(biom_16s_qc) > 0) >= 5, ] ### PCR errors  

# removes any OTUs that has less than 5 total reads across all samples
otu_table(physeq_bac_filter) <- otu_table(physeq_bac_filter)[which(rowSums(otu_table(physeq_bac_filter)) >= 5),]### PCR Errors 
otu_table(physeq_bac_filter)
tax_table(physeq_bac_filter)
sample_data(physeq_bac_filter)
physeq_bac_filter

sums_physeq_bac_filter <- data.frame(colSums(otu_table(physeq_bac_filter)))
colnames(sums_physeq_bac_filter) <- "Sample_TotalSeqs"
sums_physeq_bac_filter$Sample <- row.names(sums_physeq_bac_filter)
sums_physeq_bac_filter

#remove  samples with 0 AND MOCK
#write.csv(otu_table(physeq_bac_filter), file = "otu_filtered.csv")
#otu_table(physeq_bac_filter) <- subset(otu_table(physeq_bac_filter),
#select = -c(baGCRep1, baGCRep2, baGCRep3, baGCRep4, baGCRep5, baFCRep1, baFCRep2, baFCRep3, baFCRep4, baFCRep5 ))
#physeq_bac_filter

ggplot(sums_physeq_bac_filter, aes(x=Sample_TotalSeqs)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_TotalSeqs, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)
physeq_bac_filter


### normalizing with metagenomeseq------------------------------------- this is the package for normalizing without rarefaction
#BiocManager::install("RVAideMemoire")
library(metagenomeSeq)


# fitting into a Gaussian Model using metagenomeSeq-------------
physeq_bac_filter_norm<-physeq_bac_filter
otu_table(physeq_bac_filter_norm)
physeq_bac_normalise<-phyloseq_to_metagenomeSeq(physeq_bac_filter_norm)
p_biom<-cumNormStatFast(physeq_bac_normalise)
biom_quant<-cumNorm(physeq_bac_normalise, p=p_biom)
biom_quant
normFactors(biom_quant)
physeq_bac_normalise<-MRcounts(biom_quant, norm=T)
head(physeq_bac_normalise)
physeq_bac_normalise
#create physeq object with normalized otu table
otu_table(physeq_bac_filter_norm) <- otu_table(physeq_bac_normalise,taxa_are_rows=T)
otu_table(physeq_bac_filter_norm)

#physeq_obj_ITS_uparse_R1_mSeq <- physeq_obj_ITS_uparse_R1_clean
#otu_table(physeq_obj_ITS_uparse_R1_mSeq) <- otu_table(biom_ITS_soil, taxa_are_rows=TRUE)

physeq_bac_filter_norm
head(otu_table(physeq_bac_filter_norm))
head(tax_table(physeq_bac_filter_norm))
head(sample_data(physeq_bac_filter_norm))

write.csv(otu_table(physeq_bac_filter_norm), file = "filtered2_otus.csv")




#######################################################################################
##################           VISUALISATION             ################################
#######################################################################################
#https://joey711.github.io/phyloseq/plot_ordination-examples.html


######################   UNCONSTRANED ORDINATION  ##########################
#######################################################################################
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#unconstrained_ordinations

#PCOA

#One of the best exploratory analyses for amplicon data is unconstrained ordinations. 
#Here we will look at ordinations of our full community samples. We will use the scale_reads() function in miseqR.R to scale to the smallest library size, which is the default. 
#If you want to scale to another depth, you can do so by setting the "n" argument

theme_set(theme_bw())
# Scale reads to even depth 

physeq_bac_scale <- transform_sample_counts(physeq_bac_filter_norm, function(x) 1E6 * x/sum(x) ) 
#Keep only the most abundant fifteen phyla.
#Genus.sum = tapply(taxa_sums(physeq_bac_scale), tax_table(physeq_bac_scale)[, "Genus"], sum, na.rm=TRUE)
#top19taxa = names(sort(Genus.sum, TRUE))[1:20]
#physeq_bac_scale = prune_taxa((tax_table(physeq_bac_scale)[, "Genus"] %in% top19taxa), physeq_bac_scale)

#physeq_bac_scale.ord <- ordinate(physeq_bac_scale, "NMDS", "bray")
#p1 = plot_ordination(physeq_bac_scale, physeq_bac_scale.ord, type="taxa", color="Genus", title="taxa")
#print(p1)

# Ordinate
ordination_pcoa <- ordinate(
  physeq = physeq_bac_scale, 
  method = "PCoA", 
  distance = "bray")#,weighted=TRUE)

# Plot 
#colors <- c("Field2" = "#7FC97F", "Field4" = "#BEAED4", "Field5" = "#FDC086", "Field6" = "#386CB0", "Field8" = "#F0027F", "Field10" = "#BF5B17", "Field13" = "#E5D8BD", "Field14" = "#1B1E77", "Field15" = "#FF0000")
colors <- c("Field13" = "#7FC97F", "Field8" = "#BEAED4", "Field2" = "#FDC086", "Field4" = "#386CB0", "Field5" = "#F0027F", "Field6" = "#BF5B17", "Field10" = "#E5D8BD", "Field14" = "#1B9E11", "Field15" = "#FF0000")
#SFW <- c("Field13" = "Degraded", "Field8" = "Degraded", "Field2" = "Disturbed", "Field4" = "Disturbed", "Field5" = "Degraded", "Field6" = "Disturbed", "Field10" = "Disturbed", "Field14" = "Degraded", "Field15" = "Degraded")

shapes <- c("Field2" = "0", "Field4" = "15", "Field5" = "15", "Field6" = "15", "Field8" = "17", "Field10" = "17", "Field13" = "1", "Field14" = "1", "Field15" = "16")


p <- plot_ordination(
  physeq = physeq_bac_scale,
  ordination = ordination_pcoa,
  color = "Population",
  shape = "Soil",  
  title = "PCoA of Bacteria Communities"
) 
p = p + geom_point(aes(colour = factor(Population)))
p = p + scale_colour_manual(values = colors)
p = p + geom_point(size=7, alpha=10)
p = p + stat_ellipse(aes(group=Soilfoodweb), type="norm", alpha=4, linetype = 2, level = 0.70, show.legend = TRUE)

print(p)                     



# PERMANOVA------------------------------------------------------------------------------------------

options(scipen = 999) 
library("vegan")
library("RVAideMemoire")
set.seed(1)

# Calculate bray curtis distance matrix
physeq_bac_scale_bray <- phyloseq::distance(physeq_bac_scale, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq_field))

# Adonis test
#  a*b means a cross of a and b. This os the same as a + b + a:b
#  a:b means interactions of all terms in "a" with all terms in "b".
#  a+b means all the terms in "a" together with all the terms in "b" with duplicates removed 

factor <- model.matrix(~ Soil + Region + Soilfoodweb , data = sampledf) # do this see which factors will be compared in adonis


adonis(physeq_bac_scale_bray ~ Soil + Region + Soilfoodweb, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soil:Soilfoodweb, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soil:Region, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soilfoodweb:Region, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soil:Region:Soilfoodweb, data = sampledf, permutations=9999)

# Homogeneity of dispersion test

beta_Soil <- betadisper(physeq_bac_scale_bray, sampledf$Soil)
permutest(beta_Soil)

beta_Region <- betadisper(physeq_bac_scale_bray, sampledf$Region)
permutest(beta_Region)

beta_Soilfoodweb <- betadisper(physeq_bac_scale_bray, sampledf$Soilfoodweb)
permutest(beta_Soilfoodweb)

beta_Soil_Soilfoodweb <- betadisper(physeq_bac_scale_bray, sampledf$Soil_Soilfoodweb)
permutest(beta_Soil_Soilfoodweb)

beta_Soil_Region <- betadisper(physeq_bac_scale_bray, sampledf$Soil_Region)
permutest(beta_Soil_Region)

beta_Soilfoodweb_Region <- betadisper(physeq_bac_scale_bray, sampledf$Soilfoodweb_Region)
permutest(beta_Soilfoodweb_Region)

beta_Soil_Region_Soilfoodweb <- betadisper(physeq_bac_scale_bray, sampledf$Soil_Region_Soilfoodweb)
permutest(beta_Soil_Region_Soilfoodweb)

## CORE MICROBIOME FOR FIELD POPULATIONS---------------------------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
source("sncm.fit.R")
source("ExtractCore.R")
theme_set(theme_light())

# Remove sample with less than 10 reads 
physeq_field_filt <- prune_samples(sample_sums(physeq_field)>=10, physeq_field)

colSums(physeq_field_filt@otu_table)

# Core microbiome for field and greenhouse populations
ExtractCore(physeq_field_filt, "Population", "Increase", Group=NULL, Level=NULL) -> core_field
core_field_list<- core_field[1]

write.table(core_field_list, file = "foo.txt", sep = ",", col.names = NA,qmethod = "double") #column for otu labelled "taxa"
core_field_list <- read.table("foo.txt", header = TRUE, sep = ",", row.names = 1)

#Phyloseq with only field core

subset_core_field <- subset(otu_table(physeq_field_filt), rownames(otu_table(physeq_field_filt)) %in% (core_field_list$Taxa))
physeq_field.core <- merge_phyloseq(subset_core_field, tax_table(physeq_field_filt), sample_data(physeq_field_filt))

#Neutral model
#install_github("DanielSprockett/tyRa")
library(tyRa)

spp.out <- tyRa::fit_sncm(spp = otu_table(physeq_field.core)@.Data, pool=NULL, taxon=data.frame(tax_table(physeq_field.core)))
plot_sncm_fit(spp.out, fill = NULL, title = "Model Fit")


# Piechart of field population core microbiome----------------------------------------------------------------------------------------------------

#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(physeq_field.core, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(physeq_field.core)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)
#write.table(physeq_bacRm,"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus
pm_phylum

physeq_field.core_rel <- data.table(pm_phylum)
physeq_field.core_rel [(Abundance <= 0.000), Genus:= "Other"]

ggplot(physeq_field.core_rel, aes(x = "", y = Abundance, fill = Genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("Amycolatopsis" ="#CD9BCD", 
                               "Bradyrhizobium" ="#652926", 
                               "Candidatus Phytoplasma" = "#673770",
                               "Clostridium" ="orange" ,
                               "Pseudomonas" ="plum4" ,
                               "Rheinheimera" ="darkorchid1" ,
                               "Unclassified" ="#D9D9D9" ,
                               "Other" =	"darksalmon"
  )) +
  theme_void() + guides(fill=guide_legend(nrow=10)) 




































#remove  samples with less than 10 otus 

#otu_table(physeq_bac_filter) <- subset(otu_table(physeq_bac_filter),
#select = -c(baGCRep1, baGCRep2, baGCRep3, baGCRep4, baGCRep5, baFCRep1, baFCRep2, baFCRep3, baFCRep4, baFCRep5 ))
#physeq_bac_filter


# removes any OTUs that has less than 5 total reads across all samples
physeq_field_filt <- physeq_field
otu_table(physeq_field_filt) <- otu_table(physeq_field_filt)[which(rowSums(otu_table(physeq_field_filt)) >= 5),]### PCR Errors 

physeq_field_filt

# Core microbiome for field and greenhouse populations

ExtractCore(physeq_filt, "Origin", "Increase", Group=NULL, Level=NULL) -> core_field_greenhouse
core_field_greenhouse



physeq_filt@sam_data
colSums(physeq_field_filt@otu_table)



new_obj[4]
new_obj[6]

otu_table(physeq_field) physeq_field

FitNeutral(new_obj)


taxon <- new_obj[[7]]
spp <- t(new_obj[[5]])
colSums(spp)
sncm.fit(spp, taxon, stats = TRUE, pool = NULL)

































#******************************* GREENHOUSE MHAPLA POPULATIONS *********************************##

#RAW ABUNDANCE

#Phyla
physeq_bacRm = merge_samples(physeq_greenhouse, "Population2")
sample_data(physeq_bacRm)$ Field <- levels(sample_data(physeq_filt)$Population2)

pm_Phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate taxa at order level
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Phylum)           # Sort data frame alphabetically by Phylum
pm_Phylum

write.table(pm_Phylum, "untransformed_greenhouse_16S_phyla.txt")



## GENERA
physeq_bacRm = merge_samples(physeq_greenhouse, "Population2")
sample_data(physeq_bacRm)$ Field <- levels(sample_data(physeq_greenhouse)$Population2)

pm_Genus <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Population2)           # Sort data frame alphabetically by Genus
pm_Genus

write.table(pm_Genus, "eg.untransformed_16S_greenhouse_genera.txt")






### RELATIVE ABUNDANCE GENUS ###

#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(physeq_greenhouse, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(physeq_greenhouse)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)

pm_phylum.greenhouse <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus


Abundance.greenhouse <- data.table(pm_phylum.greenhouse)
Abundance.greenhouse [(Abundance <= 0.00)] #Abundance.greenhouse [(Abundance <= 0.00), Genus:= "Other"]

positions <- c("Greenhouse13","Greenhouse8","Greenhouse2","Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
#positions <- c("Field2","Field4","Field5","Field6","Field8","Field10","Field13","Field14","Field15")

p <- ggplot(data=Abundance.greenhouse, aes(x=Sample, y=Abundance, fill=Genus)) + theme_gray()
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("Acidibacter" = "#652926",
                               "Acidothermus" ="steelblue", 
                               "Actinophytocola" ="#C84248", 
                               "Actinospica" ="green",
                               "Amycolatopsis" ="#E5D8BD", 
                               "Bryobacter" ="#673770", 
                               "Bosea" ="#D14285", 
                               "Bradyrhizobium" ="#952926", 
                               "Brevundimonas" ="Blue",
                               "Candidatus Phytoplasma" = "#673770",
                               "Candidatus Udaeobacter" = "#AD6F3B", 
                               "Catenulispora" ="#CBD588", 
                               "Caulobacter" ="#5F7FC7", 
                               "Cellvibrio" ="orange" ,
                               "Chitinophaga" ="#DA5724",
                               "Chryseolinea" = "#508578", 
                               "Clostridium sensu stricto 1" ="#CD9BCD", 
                               "Devosia" ="tan", 
                               "Dokdonella" ="magenta1",
                               "Dongia" ="gray52" ,
                               "Duganella" ="darkorange4" ,
                               "Dyella" ="lightsalmon1" ,
                               "Flavobacterium" ="palevioletred3" ,
                               "Fluviicola" ="olivedrab1" ,
                               "Frankia" ="cyan1", 
                               "Gemmata" ="darkviolet" ,
                               "Halomonas" ="lavender", 
                               "Haliangium" ="greenyellow", 
                               "Holdemanella" ="sienna1",
                               "Hyphomicrobium" ="deepskyblue1" ,
                               "Inquilinus" ="honeydew2",
                               "Kibdelosporangium" = "red", 
                               "Kribbella" ="aquamarine4" ,
                               "Labrys" ="gold" ,
                               "Lechevalieria" ="plum3" ,
                               "Limnohabitans" ="peachpuff3", 
                               "Massilia" ="turquoise3" ,
                               "Mesorhizobium" ="blueviolet" ,
                               "Mucilaginibacter" ="green4", 
                               "Mycobacterium" ="palegreen2",
                               "Mycoplasma" ="#CCCCCC", 
                               "Niastella" ="hotpink3" ,
                               "Nocardia" ="gray38" ,
                               "Nocardioides" ="black" ,
                               "Novosphingobium" ="midnightblue" ,
                               "Ohtaekwangia" ="khaki3" ,
                               "Paenibacillus" ="yellow4" ,
                               "Paraburkholderia" ="pink",
                               "Pedomicrobium" ="cadetblue" ,
                               "Pedobacter"="purple4",
                               "Phenylobacterium" ="rosybrown3" ,
                               "Polaromonas" ="gray69" ,
                               "Pseudomonas" ="paleturquoise4",
                               "Pseudonocardia" = "limegreen" ,
                               "Rheinheimera" ="lemonchiffon3",
                               "Rhizobacter" = "deeppink1" ,
                               "Rhizobium" ="#AD6F3B",
                               "Rhodomicrobium" ="plum4" ,
                               "Rhodoplanes" ="darkolivegreen1" ,
                               "Roseiarcus" ="darkorchid1" ,
                               "SM1A02" ="rosybrown4", 
                               "Solirubrobacter" ="cornflowerblue", 
                               "Sphingobium" ="lightgoldenrodyellow", 
                               "Sphingomonas" ="springgreen1" ,
                               "Steroidobacter" ="moccasin",
                               "Streptomyces" ="lavenderblush", 
                               "Unclassified" =	"darksalmon",
                               "Variovorax" ="navyblue" ,
                               "Vibrio" ="yellow" ,
                               "Xanthomonas" ="#D9D9D9" 
  )) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=25)) + scale_x_discrete(limits = positions) + theme_gray()










### RELATIVE ABUNDANCE PHYLA ###

physeq_bacRm = merge_samples(physeq_greenhouse, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(physeq_greenhouse)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)

pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate taxa at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Phylum)           # Sort data frame alphabetically by Phylum


dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00)] #dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Phylum:= "Other"]
write.table(dat_Soil_Pre_Abs_bp,"Phylum_rel_abun_greenhouse.csv") #output the rel. abu table for further analysis in excel


positions <- c("Greenhouse13","Greenhouse8","Greenhouse2","Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
#positions <- c("Field13","Field8","Field2","Field4","Field5","Field6","Field10","Field14","Field15")

p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Phylum))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("Acidobacteriota" = "#7FC97F", 
                               "Actinobacteriota" = "#BEAED4", 
                               "Bacteroidota" = "#FDC086",
                               "Chloroflexi" = "#386CB0",
                               "Firmicutes" = "#F0027F", 
                               "Myxococcota" = "#BF5B17", 
                               "Planctomycetota" = "#666666", 
                               "Proteobacteria" = "#1B9E77", 
                               "Verrucomicrobiota" = "#FF7F00" 
  )) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) + scale_x_discrete(limits = positions) + theme_gray()






#RELATIVE ABUNDANCE BY FUNCTION


#launch_evo_palette()
#palette_box()

physeq_bacRm = merge_samples(physeq_greenhouse_func, "Population")
sample_data(physeq_bacRm)
#sample_data(physeq_bacRm)$Population <- levels(sample_data(physeq_field_func)$Population)
#physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
#View(physeq_bacRm)

pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Phylum


dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00)] #dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Phylum:= "Other"]

positions <- c("Greenhouse13","Greenhouse8","Greenhouse2","Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
#positions <- c("Field13","Field8","Field2","Field4","Field5","Field6","Field10","Field14","Field15")


p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Genus))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) + 
  scale_fill_manual(values = c("Animal_pathogenic" ="deeppink1",
                               "Anti_bacteria" = "#652926",
                               "Anti_fungi" ="steelblue",
                               "Antibiotic resistant" ="#E5D8BD",
                               "Enhanced_nematode_parasitism" ="#C84248", 
                               "Iron_reducing" ="#5F7FC7", 
                               "Nematicidal" ="Blue",
                               "Nitrogen_fixer" ="green",
                               "Other" = "#673770", 
                               "Plant_growth_promoter" ="#DA5724",
                               "Plant_pathogenic" = "#AD6F3B",
                               "Polysaccharide_degrader" ="#673770", 
                               "Root_knot_nematode_associated" ="#CBD588",
                               "Soybean_cyst_associated" ="#D14285", 
                               "Suppressive_soils" ="orange"
  )) +
  theme(legend.position="right") + scale_x_discrete(limits = positions) + theme_gray() + guides(fill=guide_legend(nrow=15))


















#Splitting populations to identify the presence/absence of genera-----------------------------------------------------------------
Abundance.greenhouse.split <- split(Abundance.greenhouse, Abundance.greenhouse$Sample)

#Population 13
Abundance.greenhouse.split.13 <- Abundance.greenhouse.split$Greenhouse13
Abundance.greenhouse.split.13.genera <- unique(Abundance.greenhouse.split.13$Genus)

#Population 8
Abundance.greenhouse.split.8 <- Abundance.greenhouse.split$Greenhouse8
Abundance.greenhouse.split.8.genera <- unique(Abundance.greenhouse.split.8$Genus)

#Population 2
Abundance.greenhouse.split.2 <- Abundance.greenhouse.split$Greenhouse2
Abundance.greenhouse.split.2.genera <- unique(Abundance.greenhouse.split.2$Genus)

#Population 4
Abundance.greenhouse.split.4 <- Abundance.greenhouse.split$Greenhouse4
Abundance.greenhouse.split.4.genera <- unique(Abundance.greenhouse.split.4$Genus)

#Population 5
Abundance.greenhouse.split.5 <- Abundance.greenhouse.split$Greenhouse5
Abundance.greenhouse.split.5.genera <- unique(Abundance.greenhouse.split.5$Genus)

#Population 6
Abundance.greenhouse.split.6 <- Abundance.greenhouse.split$Greenhouse6
Abundance.greenhouse.split.6.genera <- unique(Abundance.greenhouse.split.6$Genus)

#Population 10
Abundance.greenhouse.split.10 <- Abundance.greenhouse.split$Greenhouse10
Abundance.greenhouse.split.10.genera <- unique(Abundance.greenhouse.split.10$Genus)

#Population 14
Abundance.greenhouse.split.14 <- Abundance.greenhouse.split$Greenhouse14
Abundance.greenhouse.split.14.genera <- unique(Abundance.greenhouse.split.14$Genus)

#Population 15
Abundance.greenhouse.split.15 <- Abundance.greenhouse.split$Greenhouse15
Abundance.greenhouse.split.15.genera <- unique(Abundance.greenhouse.split.15$Genus)
Abundance.greenhouse.split.15.genera <- unique(Abundance.greenhouse.split.15$Genus)

#All genera in greenhouse populations
greenhouse.genera <- unique(Abundance.greenhouse$Genus)

# Finding maximum length
max_ln1 <- max(c(length(Abundance.greenhouse.split.13.genera), length(Abundance.greenhouse.split.8.genera)))
max_ln2 <- max(c(length(Abundance.greenhouse.split.2.genera), length(Abundance.greenhouse.split.4.genera)))
max_ln3 <- max(c(length(Abundance.greenhouse.split.5.genera), length(Abundance.greenhouse.split.6.genera)))
max_ln4 <- max(c(length(Abundance.greenhouse.split.10.genera), length(Abundance.greenhouse.split.14.genera)))
max_ln5 <- max(c(length(Abundance.greenhouse.split.15.genera), length(greenhouse.genera)))

max_ln<-max(max_ln1,max_ln2,max_ln3,max_ln4,max_ln5)

# Merge all populations
Pop.greenhouse.13 = c(Abundance.greenhouse.split.13.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.13.genera)))
Pop.greenhouse.8 = c(Abundance.greenhouse.split.8.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.8.genera)))
Pop.greenhouse.2 = c(Abundance.greenhouse.split.2.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.2.genera)))
Pop.greenhouse.4 = c(Abundance.greenhouse.split.4.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.4.genera)))
Pop.greenhouse.5 = c(Abundance.greenhouse.split.5.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.5.genera)))
Pop.greenhouse.6 = c(Abundance.greenhouse.split.6.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.6.genera)))
Pop.greenhouse.10 = c(Abundance.greenhouse.split.10.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.10.genera)))
Pop.greenhouse.14 = c(Abundance.greenhouse.split.14.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.14.genera)))
Pop.greenhouse.15 = c(Abundance.greenhouse.split.15.genera,rep(NA, max_ln - length(Abundance.greenhouse.split.15.genera)))
All.greenhouse.genera = c(greenhouse.genera,rep(NA, max_ln - length(greenhouse.genera)))


Combined.greenhouse.genera <- cbind(All.greenhouse.genera, Pop.greenhouse.13, Pop.greenhouse.8, Pop.greenhouse.2, Pop.greenhouse.4, Pop.greenhouse.5, Pop.greenhouse.6, Pop.greenhouse.10, Pop.greenhouse.14, Pop.greenhouse.15)

write.csv(Combined.greenhouse.genera, file="All.greenhouse.genus2.csv")



###    ALPHA DIVERSITY ###

#load libraries
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(vegan)
library(ggplot2)

#plot alpha diversity
ps1 <- prune_taxa(taxa_sums(Mollicutes.greenhouse) > 0, Mollicutes.greenhouse) ## use unfiltered phyloseq object for alpha diversity
alpha_diversity <- alpha(ps1, index = "diversity_shannon")
head(alpha_diversity)


#Prepare data for visualisation
ps1.meta <- meta(ps1)
head(ps1.meta)

#Add the diversity table to metadata
ps1.meta$diversity_observed  <- alpha_diversity$observed 
ps1.meta$diversity_shannon <- alpha_diversity$diversity_shannon

#Observed diversity/Richness
ggplot(ps1.meta, aes(x = Population, y = diversity_observed)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.
kruskal.test(diversity_observed ~ Population, data=ps1.meta) #p-value = 0.3328


#Shannon diversity
ggplot(ps1.meta, aes(x = Population, y = diversity_shannon)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.
kruskal.test(diversity_shannon ~ Population, data=ps1.meta) #p-value = 0.3953


###    ALPHA DIVERSITY MOLLICUTE   ###

Mollicutes.greenhouse <- subset_taxa(physeq_greenhouse, Class=="Mollicutes")

#plot alpha diversity
ps1 <- prune_taxa(taxa_sums(Mollicutes.greenhouse) > 0, Mollicutes.greenhouse) ## use unfiltered phyloseq object for alpha diversity
alpha_diversity <- alpha(ps1, index = "all")
head(alpha_diversity)


#Prepare data for visualisation
ps1.meta <- meta(ps1)
head(ps1.meta)

#Add the diversity table to metadata
ps1.meta$diversity_observed  <- alpha_diversity$diversity_observed 
ps1.meta$diversity_shannon <- alpha_diversity$diversity_shannon

#Observed diversity/Richness
ggplot(ps1.meta, aes(x = Population, y = diversity_observed)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.
kruskal.test(diversity_observed ~ Population, data=ps1.meta) #p-value = 0.3328


#Shannon diversity
ggplot(ps1.meta, aes(x = Population, y = diversity_shannon)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.
kruskal.test(diversity_shannon ~ Population, data=ps1.meta) #p-value = 0.3953












###   BETA DIVERSITY   ####

# any sample with less than 5 reads for a particular otu will be placed to 0
physeq_bac_filter <- physeq_greenhouse
#otu_table(physeq_bac_filter)[otu_table(physeq_bac_filter) <= 5] < 0- 0 ### tag switching
#otu_table(biom_16s_qc) <- otu_table(biom_16s_qc)[rowSums(otu_table(biom_16s_qc) > 0) >= 5, ] ### PCR errors  

# removes any OTUs that has less than 5 total reads across all samples
otu_table(physeq_bac_filter) <- otu_table(physeq_bac_filter)[which(rowSums(otu_table(physeq_bac_filter)) >= 5),]### PCR Errors 
otu_table(physeq_bac_filter)
tax_table(physeq_bac_filter)
sample_data(physeq_bac_filter)
physeq_bac_filter

sums_physeq_bac_filter <- data.frame(colSums(otu_table(physeq_bac_filter)))
colnames(sums_physeq_bac_filter) <- "Sample_TotalSeqs"
sums_physeq_bac_filter$Sample <- row.names(sums_physeq_bac_filter)
sums_physeq_bac_filter

#remove  samples with 0 AND MOCK
#write.csv(otu_table(physeq_bac_filter), file = "otu_filtered.csv")
#otu_table(physeq_bac_filter) <- subset(otu_table(physeq_bac_filter),
#select = -c(baGCRep1, baGCRep2, baGCRep3, baGCRep4, baGCRep5, baFCRep1, baFCRep2, baFCRep3, baFCRep4, baFCRep5 ))
#physeq_bac_filter

ggplot(sums_physeq_bac_filter, aes(x=Sample_TotalSeqs)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_TotalSeqs, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)
physeq_bac_filter


### normalizing with metagenomeseq------------------------------------- this is the package for normalizing without rarefaction
#BiocManager::install("RVAideMemoire")
library(metagenomeSeq)


# fitting into a Gaussian Model using metagenomeSeq-------------
physeq_bac_filter_norm<-physeq_bac_filter
otu_table(physeq_bac_filter_norm)
physeq_bac_normalise<-phyloseq_to_metagenomeSeq(physeq_bac_filter_norm)
p_biom<-cumNormStatFast(physeq_bac_normalise)
biom_quant<-cumNorm(physeq_bac_normalise, p=p_biom)
biom_quant
normFactors(biom_quant)
physeq_bac_normalise<-MRcounts(biom_quant, norm=T)
head(physeq_bac_normalise)
physeq_bac_normalise
#create physeq object with normalized otu table
otu_table(physeq_bac_filter_norm) <- otu_table(physeq_bac_normalise,taxa_are_rows=T)
otu_table(physeq_bac_filter_norm)

#physeq_obj_ITS_uparse_R1_mSeq <- physeq_obj_ITS_uparse_R1_clean
#otu_table(physeq_obj_ITS_uparse_R1_mSeq) <- otu_table(biom_ITS_soil, taxa_are_rows=TRUE)

physeq_bac_filter_norm
head(otu_table(physeq_bac_filter_norm))
head(tax_table(physeq_bac_filter_norm))
head(sample_data(physeq_bac_filter_norm))

write.csv(otu_table(physeq_bac_filter_norm), file = "filtered2_otus.csv")




#######################################################################################
##################           VISUALISATION             ################################
#######################################################################################
#https://joey711.github.io/phyloseq/plot_ordination-examples.html


######################   UNCONSTRANED ORDINATION  ##########################
#######################################################################################
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#unconstrained_ordinations

#PCOA

#One of the best exploratory analyses for amplicon data is unconstrained ordinations. 
#Here we will look at ordinations of our full community samples. We will use the scale_reads() function in miseqR.R to scale to the smallest library size, which is the default. 
#If you want to scale to another depth, you can do so by setting the "n" argument

theme_set(theme_bw())
# Scale reads to even depth 

physeq_bac_scale <- transform_sample_counts(physeq_bac_filter_norm, function(x) 1E6 * x/sum(x) ) 
#Keep only the most abundant fifteen phyla.
#Genus.sum = tapply(taxa_sums(physeq_bac_scale), tax_table(physeq_bac_scale)[, "Genus"], sum, na.rm=TRUE)
#top19taxa = names(sort(Genus.sum, TRUE))[1:20]
#physeq_bac_scale = prune_taxa((tax_table(physeq_bac_scale)[, "Genus"] %in% top19taxa), physeq_bac_scale)

#physeq_bac_scale.ord <- ordinate(physeq_bac_scale, "NMDS", "bray")
#p1 = plot_ordination(physeq_bac_scale, physeq_bac_scale.ord, type="taxa", color="Genus", title="taxa")
#print(p1)

# Ordinate
ordination_pcoa <- ordinate(
  physeq = physeq_bac_scale, 
  method = "PCoA", 
  distance = "bray")#,weighted=TRUE)

# Plot 
colors <- c("Greenhouse13" = "#7FC97F", "Greenhouse8" = "#BEAED4", "Greenhouse2" = "#FDC086", "Greenhouse4" = "#386CB0", "Greenhouse5" = "#F0027F", "Greenhouse6" = "#BF5B17", "Greenhouse10" = "#E5D8BD", "Greenhouse14" = "#1B9E11", "Greenhouse15" = "#FF0000")

shapes <- c("Greenhouse2" = "0", "Greenhouse4" = "15", "Greenhouse5" = "15", "Greenhouse6" = "15", "Greenhouse8" = "17", "Greenhouse10" = "17", "Greenhouse13" = "1", "Greenhouse14" = "1", "Greenhouse15" = "16")


p <- plot_ordination(
  physeq = physeq_bac_scale,
  ordination = ordination_pcoa,
  color = "Population",
  shape = "Soil",  title = "PCoA of Bacteria Communities"
) 
p = p + geom_point(aes(colour = factor(Population)))
p = p + scale_colour_manual(values = colors)
p = p + geom_point(size=7, alpha=10)
p = p + stat_ellipse(aes(group=Soilfoodweb), type="norm", alpha=4, linetype = 2, level = 0.70, show.legend = TRUE)

print(p)                     




# PERMANOVA------------------------------------------------------------------------------------------

options(scipen = 999) 
library("vegan")
library("RVAideMemoire")
set.seed(1)

# Calculate bray curtis distance matrix
physeq_bac_scale_bray <- phyloseq::distance(physeq_bac_scale, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq_field))

# Adonis test
#  a*b means a cross of a and b. This os the same as a + b + a:b
#  a:b means interactions of all terms in "a" with all terms in "b".
#  a+b means all the terms in "a" together with all the terms in "b" with duplicates removed 

factor <- model.matrix(~ Soil + Region + Soilfoodweb , data = sampledf) # do this see which factors will be compared in adonis


adonis(physeq_bac_scale_bray ~ Soil + Region + Soilfoodweb, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soil:Soilfoodweb, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soil:Region, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soilfoodweb:Region, data = sampledf, permutations=9999)
adonis(physeq_bac_scale_bray ~ Soil:Region:Soilfoodweb, data = sampledf, permutations=9999)




# Homogeneity of dispersion test

beta_Soil <- betadisper(physeq_bac_scale_bray, sampledf$Soil)
permutest(beta_Soil)

beta_Region <- betadisper(physeq_bac_scale_bray, sampledf$Region)
permutest(beta_Region)

beta_Soilfoodweb <- betadisper(physeq_bac_scale_bray, sampledf$Soilfoodweb)
permutest(beta_Soilfoodweb)

beta_Soil_Soilfoodweb <- betadisper(physeq_bac_scale_bray, sampledf$Soil_Soilfoodweb)
permutest(beta_Soil_Soilfoodweb)

beta_Soil_Region <- betadisper(physeq_bac_scale_bray, sampledf$Soil_Region)
permutest(beta_Soil_Region)

beta_Soilfoodweb_Region <- betadisper(physeq_bac_scale_bray, sampledf$Soilfoodweb_Region)
permutest(beta_Soilfoodweb_Region)

beta_Soil_Region_Soilfoodweb <- betadisper(physeq_bac_scale_bray, sampledf$Soil_Region_Soilfoodweb)
permutest(beta_Soil_Region_Soilfoodweb)

















## CORE MICROBIOME FOR FIELD POPULATIONS---------------------------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(vegan)
library(tidyr)
library(dplyr)
library(ggplot2)
source("sncm.fit.R")
source("ExtractCore.R")
theme_set(theme_light())

# Remove sample with less than 10 reads 
physeq_greenhouse_filt <- prune_samples(sample_sums(physeq_greenhouse)>=10, physeq_greenhouse)

colSums(physeq_greenhouse_filt@otu_table)

# Core microbiome for field and greenhouse populations
ExtractCore(physeq_greenhouse_filt, "Population", "Increase", Group=NULL, Level=NULL) -> core_field
core_field_list<- core_field[1]

write.table(core_field_list, file = "foo1.txt", sep = ",", col.names = NA,qmethod = "double") #column for otu labelled "taxa"
core_field_list <- read.table("foo1.txt", header = TRUE, sep = ",", row.names = 1)

#Phyloseq with only field core

subset_core_field <- subset(otu_table(physeq_greenhouse_filt), rownames(otu_table(physeq_greenhouse_filt)) %in% (core_field_list$Taxa))
physeq_greenhouse.core <- merge_phyloseq(subset_core_field, tax_table(physeq_greenhouse_filt), sample_data(physeq_greenhouse_filt))

#Neutral model
#install_github("DanielSprockett/tyRa")
library(tyRa)

spp.out <- tyRa::fit_sncm(spp = otu_table(physeq_greenhouse.core)@.Data, pool=NULL, taxon=data.frame(tax_table(physeq_greenhouse.core)))
plot_sncm_fit(spp.out, fill = NULL, title = "Model Fit")


# Piechart of field population core microbiome----------------------------------------------------------------------------------------------------

#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(physeq_greenhouse.core, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(physeq_greenhouse.core)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)
#write.table(physeq_bacRm,"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus


physeq_greenhouse.core_rel <- data.table(pm_phylum)
physeq_greenhouse.core_rel [(Abundance <= 0.000), Genus:= "Other"]

ggplot(physeq_greenhouse.core_rel, aes(x = "", y = Abundance, fill = Genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("Amycolatopsis" ="#CD9BCD", 
                               "Bradyrhizobium" ="#652926", 
                               "Candidatus Phytoplasma" = "#673770",
                               "Cellvibrio" ="#CBD588", 
                               "Clostridium" ="orange" ,
                               "Pseudomonas" ="plum4" ,
                               "Pseudonocardia" ="darkolivegreen1" ,
                               "Rheinheimera" ="darkorchid1" ,
                               "Unclassified" ="#D9D9D9" ,
                               "Other" =	"darksalmon"
  )) +
  theme_void() + guides(fill=guide_legend(nrow=10)) 

























## CORE MICROBIOME FOR GREENHOUSE POPULATIONS---------------------------------------------------------------------------------------------

library(tidyverse)
library(reshape2)
library(vegan)
library(tidyr)
library(dplyr)
theme_set(theme_light())


nReads=10000                                                            # input dataset needs to be rarified and the rarifaction depth included 

#write.table(otu_table(physeq_filt),file="otu.txt",sep='\t', col.names=
#TRUE, row.names = TRUE)

#otu<- read.table(file="otu.txt",header=TRUE)

otu <- as.data.frame(otu_table(physeq_greenhouse)) 
head(otu)
map <- as.data.frame(sample_data(physeq_greenhouse))
head(map)

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance  
occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance 

PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sample_ID, abun, -otu) %>%
  left_join(map, by = 'sample_ID') %>%
  group_by(otu, Population) %>%
  summarise(plot_freq=sum(abun>0)/length(abun),        # frequency of detection between time points
            coreSite=ifelse(plot_freq == 1, 1, 0), # 1 only if occupancy 1 with specific genotype, 0 if not
            detect=ifelse(plot_freq > 0, 1, 0)) %>%    # 1 if detected and 0 if not detected with specific genotype
  group_by(otu) %>%
  summarise(sumF=sum(plot_freq),
            sumG=sum(coreSite),
            nS=length(Population)*2,
            Index=(sumF+sumG)/nS) # calculating weighting Index based on number of time points detected and 

otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,
            rank=Index) %>%
  arrange(desc(rank))

BCaddition <- NULL

otu_start=otu_ranked$otu[1]
start_matrix <- as.matrix(otu[otu_start,])
#start_matrix <- t(start_matrix)
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)
names(df_s)[2] <- 1 
BCaddition <- rbind(BCaddition,df_s)

# calculating BC dissimilarity based on additon of ranked OTUs from 2nd to 2067th. Can be set to the entire length of OTUs in the dataset, however it might take some Population if more than 5000 OTUs are included.
for(i in 2:675){                              
  otu_add=otu_ranked$otu[i]                       
  add_matrix <- as.matrix(otu[otu_add,])
  add_matrix <- (add_matrix)
  start_matrix <- rbind(start_matrix, add_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)
  names(df_a)[2] <- i 
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
}
# calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs) 
x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))   
x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
df_full <- data.frame(x_names,x)
names(df_full)[2] <- length(rownames(otu))
BCfull <- left_join(BCaddition,df_full, by='x_names')

rownames(BCfull) <- BCfull$x_names
temp_BC <- BCfull
temp_BC$x_names <- NULL
temp_BC_matrix <- as.matrix(temp_BC)

BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
BC_ranked <- left_join(BC_ranked, increaseDF)
BC_ranked <- BC_ranked[-nrow(BC_ranked),]

#Creating thresholds for core inclusion 

#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
  return(left - right)
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

elbow <- which.max(BC_ranked$fo_diffs)

#B) Final increase in BC similarity of equal or greater then 2% 
lastCall <- last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))

#Creating plot of Bray-Curtis similarity
ggplot(BC_ranked[1:3000,], aes(x=factor(BC_ranked$rank[1:3000], levels=BC_ranked$rank[1:3000]))) +
  geom_point(aes(y=proportionBC)) +
  theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
  geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])), lty=3, col='blue', cex=.5) +
  labs(x='ranked OTUs',y='Bray-Curtis similarity') +
  annotate(geom="text", x=elbow+14, y=.1, label=paste("Elbow method"," (",elbow,")", sep=''), color="red")+    
  annotate(geom="text", x=last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))+3, y=.5, label=paste("Last 2% increase (",last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)])),")",sep=''), color="blue")

#Creating occupancy abundance plot
occ_abun$fill <- 'no'
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.02)]))]] <- 'core'

#Fitting neutral model (Burns et al., 2016 (ISME J) - functions are in the sncm.fit.R)
source("sncm.fit.R")
spp=t(otu)
taxon=as.vector(rownames(otu))

#Models for the whole community
#obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
#sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
#sta.np.16S <- sta.np

#above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness
#below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness

#ap = obs.np$freq > (obs.np$pred.upr)
#bp = obs.np$freq < (obs.np$pred.lwr)

#ggplot() +
# geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2)+
#geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='blue', size=1.8) +
#geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +
#geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25)+
#geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25)+
#labs(x="log10(mean relative abundance)", y="Occupancy")

#Creating a plot of core taxa occupancy by Population
core <- occ_abun$otu[occ_abun$fill == 'core']

otu_relabun <- decostand(otu, method="total", MARGIN=2)

plotDF <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
  gather(sample_ID, relabun, -otu) %>%
  left_join(map, by = 'sample_ID') %>%
  left_join(otu_ranked, bu='otu') %>%
  filter(otu %in% core) %>% 
  group_by(otu, Population) %>%
  summarise(Population_freq=sum(relabun>0)/length(relabun),        
            corePopulation=ifelse(Population_freq == 1, 1, 0),      
            detect=ifelse(Population_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels=otu_ranked$otu[1:675])

ggplot(plotDF,aes(x=otu, Population_freq,fill=factor(Population))) +    
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(plotDF$otu))) +
  theme(axis.text = element_text(size=6)) +
  labs(x='Ranked OTUs', y='Occupancy by Population')

# Creating a new phyloseq with core taxa for downstream analysis

Subset.core <- subset(otu_table(physeq_filt), rownames(otu_table(physeq_filt)) %in% core)
physeq_greenhouse.core <- merge_phyloseq(Subset.core, tax_table(physeq_greenhouse), sample_data(physeq_greenhouse))


# CORE FIELD POPULATIONS PIECHART----------------------------------------------------------------------------------------------------

#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(physeq_greenhouse.core, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(physeq_greenhouse.core)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)
#write.table(physeq_bacRm,"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Genus)           # Sort data frame alphabetically by Genus
pm_phylum

dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.01), Genus:= "Other"]

#manually ordering levels
dat_Soil_Pre_Abs_bp1 <- dat_Soil_Pre_Abs_bp
dat_Soil_Pre_Abs_bp1
#positions <- c("Field1","Field2","Field3","Field4","Field5","Field6","Field7","Field8","Field9","Field10","Field11","Field12","Field13","Field14","Field15")
#positions <- c("Natural1","Natural2","Natural6","Natural7","Natural8","Natural9","Natural10","Natural11","Natural12","Natural13","Natural15")
#positions <- c("Field4","F5","F6","F10","F14","F15","F2","F8","F13","F1","F3","F7","F9","F11","F12")
#positions <- c("N6","N10","N15","N2","N8","N13","N1","N7","N9","N11","N12")

ggplot(dat_Soil_Pre_Abs_bp, aes(x = "", y = Abundance, fill = Genus)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" , "#E31A1C",
                               "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" ,"#FBB4AE", "#B3CDE3" ,"#CCEBC5", "#DECBE4", "#FED9A6",
                               "#FFFFCC" ,"#E5D8BD", "#FDDAEC", "#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
                               "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999" ,"#66C2A5",
                               "#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072",
                               "#80B1D3", "#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F","#33A02C", "#FB9A99")) +
  theme_void() 






#Mollicute










































































#scale_color_manual(values = c("#a65628", "red", "#ffae19",
"#4daf4a", "#1919ff", "darkorchid3", "magenta")

#p <- plot_ordination(physeq_bac_scale, ordu, color="Field", shape="Landscape")
#p <- p + geom_point(size=7, alpha=.7)
#p <- p + scale_colour_brewer(type="qual", palette="Set1")
#p <- p + ggtitle("MDS/PCoA on unweighted-UniFrac distance")

#print(p)


p <- ordinate(physeq_bac_scale, "PCoA", "unifrac", weighted=TRUE)
p <- plot_ordination(physeq_bac_scale, ordu, color="Field", shape="Landscape")
p <- p + geom_point(size=7, alpha=.7)
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + ggtitle("MDS/PCoA on unweighted-UniFrac distance")

print(p)





















 ###################  DATA SUMMARY   #########################
        ##################################################################################
#setting working directory
#wdir <- "C:/Users/ISAAC/Google Drive/WORK/phd/5 MANUSCRIPTS/3 Soil Microbiome/Analysis/2018NRKNPV"
output <-(paste( wdir,"/Output/", sep = ""))
source(paste( wdir, "/Alpha_div_function.r", sep = "" ))
source("miseqR.R")


#BiocManager::install("microbiome")
library(phyloseq)
library(ggplot2)## for plotting
library(magrittr)
library(ggpubr)##for combining the plots
library(vegan)##for community ecology based codes
library(limma)
library(edgeR)
library(Hmisc)
library(igraph)
library(labdsv)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(microbiome) #to summarize phyloseq
library(knitr)
library(fantaxtic) #used for creating top taxa relative abundance
library(data.table) ## Creating data table for relative abundance

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(biom_16S_uparse_filt))
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# SUMMARIze FUNGI DaTA 
summarize_phyloseq(biom_16S_uparse)

#Generate a table of selected diversity indicators
alpha_diversity <- alpha(biom_16S_uparse_filt, index = "all")
head(alpha_diversity)
#Richness
richness <- richness(biom_16S_uparse)
head(richness)
#Dominance  # Absolute abundances for the single most abundant taxa in each sample
dominance <- dominance(biom_16S_uparse_filt, index = "all")
head(dominance)
#Rarity and low abundance - The rarity indices quantify the concentration of rare or low abundance taxa. Various rarity indices are available (see the function help for a list of options).
rarity <- rarity(biom_16S_uparse_filt, index = "all")
kable(head(rarity))
#Coverage - The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).
coverage <- coverage(biom_16S_uparse_filt, threshold = 0.5)
head(coverage)
#Core abundance - The core_abundance function refers to the relative proportion of the core species. Non-core abundance provides the complement (1-x; see rare_abundance).
core_abundance <- core_abundance(biom_16S_uparse_filt, detection = .1/100, prevalence = 50/100)
head(core_abundance)
#Gini index - Gini index is a common measure for inequality in economical income. The inverse gini index (1/x) can also be used as a community diversity measure.
Gini_index <- inequality(biom_16S_uparse_filt)
#Evenness
Evenness <- evenness(biom_16S_uparse_filt, "all")
head(Evenness)


library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))



                  ####################### RELATIVE ABUNDANCE #############################
          ####################################################################################

#Relative abundance otu table to 100%
physeq_bacRm = merge_samples(biom_16S_uparse_filt, "Soil")
sample_data(physeq_bacRm)$ Soil <- levels(sample_data(biom_16S_uparse_filt)$Soil)
#physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
#View(physeq_bacRm)
write.table(sample_data(physeq_bacRm),"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Phylum)           # Sort data frame alphabetically by class
pm_phylum

dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.02), Phylum:= "Other"]
  
#pm_phylum <- get_top_taxa(dat_Soil_Pre_Abs_bp, 20, relative = TRUE, other_label = "Other")
#Plot overall abundance in ggplot 
p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Phylum)) + theme_gray()
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C",
                                "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" ,"#FBB4AE", "#B3CDE3" ,"#CCEBC5", "#DECBE4", "#FED9A6",
                               "#FFFFCC" ,"#E5D8BD", "#FDDAEC", "#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
                               "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999" ,"#66C2A5",
                               "#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072",
                               "#80B1D3", "#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F")) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=21)) 







#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(biom_16S_uparse_filt, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(biom_16S_uparse_filt)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)
#write.table(physeq_bacRm,"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_phylum <- physeq_bacRm  %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Phylum)           # Sort data frame alphabetically by Phylum
pm_phylum

dat_Soil_Pre_Abs_bp <- data.table(pm_phylum)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Phylum:= "Other"]

#manually ordering levels
dat_Soil_Pre_Abs_bp1 <- dat_Soil_Pre_Abs_bp
dat_Soil_Pre_Abs_bp1
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13","Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")

#pm_phylum <- get_top_taxa(dat_Soil_Pre_Abs_bp, 20, relative = TRUE, other_label = "Other")
#Plot overall abundance in ggplot 
p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Phylum))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C",
                               "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" ,"#FBB4AE", "#B3CDE3" ,"#CCEBC5", "#DECBE4", "#FED9A6",
                               "#FFFFCC" ,"#E5D8BD", "#FDDAEC", "#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
                               "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999" ,"#66C2A5",
                               "#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072",
                               "#80B1D3", "#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F")) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=21)) + scale_x_discrete(limits = positions)

#+ scale_x_discrete(limits = positions)+
#scale_x_discrete(limits= c("Field1","Field2","Field3","Field4","Field5","Field6","Field7","Field8","Field9","Field10","Field11","Field12","Field13","Field14","Field15","Natural1","Natural2","Natural3","Natural4","Natural5","Natural6","Natural7","Natural8","Natural9","Natural10","Natural11","Natural12","Natural13","Natural14","Natural15")
                  


                     #############################  ALPHA DIVERSITY ##################################
                 #######################################################################################

#load libraries
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

#plot alpha diversity
ps1 <- prune_taxa(taxa_sums(physeq_field) > 0, physeq_field) ## use unfiltered phyloseq object for alpha diversity
alpha_diversity <- alpha(ps1, index = "all")
head(alpha_diversity)


#Prepare data for visualisation
ps1.meta <- meta(ps1)
head(ps1.meta)

#Add the diversity table to metadata
ps1.meta$diversity_observed  <- alpha_diversity$observed 
ps1.meta$diversity_shannon <- alpha_diversity$diversity_shannon
#ps1.meta$diversity_inverse_simpson <- alpha_diversity$diversity_inverse_simpson

#Let's say we want to compare differences in Shannon index between Soil_Pre_Abs of the study subjects.
# create a list of pairwise comparisons
Nutrient <- levels(ps1.meta$SoilMhapla) # get the variables
# make a pairwise list that we want to compare.
Nutrient.pairs <- combn(seq_along(Nutrient), 2, simplify = FALSE, FUN = function(i)Nutrient[i])
#print(Nutrient.pairs)

#Create 1x1 plot environment so that we can see all 2 metrics at once. 
par(mfrow = c(1, 2))

#Violin plot
p1 <- ggboxplot(ps1.meta, x = "Field", y = "diversity_observed",
                add = "boxplot", fill = "Field") 
print(p1)

#Statistics
#p1 <- p1 + stat_compare_means(comparisons = Nutrient.pairs) 
#print(p1)
kruskal.test(diversity_observed ~ Field, data=ps1.meta) #significant 0.00000004252
pairwise.wilcox.test(ps1.meta$diversity_observed, ps1.meta$Field, p.adjust.method="fdr")

#Violin plot
p2 <- ggboxplot(ps1.meta, x = "Field", y = "diversity_shannon",
                add = "boxplot", fill = "Field") 
print(p2)

#Statistics
kruskal.test(diversity_shannon ~ Soil_Pre_Abs, data=ps1.meta) #significant 0.000004937
pairwise.wilcox.test(ps1.meta$diversity_shannon, ps1.meta$Soil_Pre_Abs, p.adjust.method="fdr")


########################For Fields######################################

#Violin plot
p3 <- ggboxplot(ps1.meta, x = "Population", y = "diversity_observed",
                add = "boxplot", fill = "Population")+ theme_bw()
print(p3)

#Statistics
#p1 <- p1 + stat_compare_means(comparisons = Nutrient.pairs) 
#print(p1)
kruskal.test(diversity_observed ~ Population, data=ps1.meta) #significant 0.00000004252
pairwise.wilcox.test(ps1.meta$diversity_observed, ps1.meta$Populationd, p.adjust.method="fdr")




#Violin plot
p4 <- ggboxplot(ps1.meta, x = "Field", y = "diversity_shannon",
                add = "boxplot", fill = "Field") 
print(p4)

#Statistics
#p1 <- p1 + stat_compare_means(comparisons = Nutrient.pairs) 
#print(p1)
kruskal.test(diversity_shannon ~ Field, data=ps1.meta) #significant 0.00000004252
pairwise.wilcox.test(ps1.meta$diversity_shannon, ps1.meta$Field, p.adjust.method="fdr")



require(gridExtra)
grid.arrange(p1, p2, ncol=2)

###More visualisation analysis https://rpkgs.datanovia.com/ggpubr/index.html

positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")
positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")

library(vegan)
library(ggplot2)
ggplot(ps1.meta, aes(x = Population, y = diversity_observed)) + 
  geom_boxplot() + scale_x_discrete(limits = positions) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw()
# We can use analysis of variance (ANOVA) to tell if at least one of the diversity means is different from the rest.

kruskal.test(diversity_observed ~ Population, data=ps1.meta) #p-value = 0.3176
#pairwise.wilcox.test(ps1.meta$diversity_observed, ps1.meta$Field, p.adjust.method="fdr")
#Or
#wilcox https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html



                      #########################  BETA DIVERSITY ####################################
                  #######################################################################################

# Filtering otus before beta diversity analysis----- meant to remove reads that mya have appeared due to tag switching
# filtering otus--------------------------------------------------------- meant to remove reads that mya have appeared due to tag switching



# any sample with less than 5 reads for a particular otu will be placed to 0
physeq_bac_filter <- physeq_bacRm
#otu_table(physeq_bac_filter)[otu_table(physeq_bac_filter) <= 5] < 0- 0 ### tag switching
#otu_table(biom_16s_qc) <- otu_table(biom_16s_qc)[rowSums(otu_table(biom_16s_qc) > 0) >= 5, ] ### PCR errors  

# removes any OTUs that has less than 5 total reads across all samples
otu_table(physeq_bac_filter) <- otu_table(physeq_bac_filter)[which(rowSums(otu_table(physeq_bac_filter)) >= 5),]### PCR Errors 
otu_table(physeq_bac_filter)
tax_table(physeq_bac_filter)
sample_data(physeq_bac_filter)
physeq_bac_filter

sums_physeq_bac_filter <- data.frame(colSums(otu_table(physeq_bac_filter)))
colnames(sums_physeq_bac_filter) <- "Sample_TotalSeqs"
sums_physeq_bac_filter$Sample <- row.names(sums_physeq_bac_filter)
sums_physeq_bac_filter

#remove  samples with 0 AND MOCK
#write.csv(otu_table(physeq_bac_filter), file = "otu_filtered.csv")
#otu_table(physeq_bac_filter) <- subset(otu_table(physeq_bac_filter),
#select = -c(baGCRep1, baGCRep2, baGCRep3, baGCRep4, baGCRep5, baFCRep1, baFCRep2, baFCRep3, baFCRep4, baFCRep5 ))
#physeq_bac_filter

ggplot(sums_physeq_bac_filter, aes(x=Sample_TotalSeqs)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_TotalSeqs, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)
physeq_bac_filter


### normalizing with metagenomeseq------------------------------------- this is the package for normalizing without rarefaction
#BiocManager::install("RVAideMemoire")
library(metagenomeSeq)


# fitting into a Gaussian Model using metagenomeSeq-------------
physeq_bac_filter_norm<-physeq_bac_filter
otu_table(physeq_bac_filter_norm)
physeq_bac_normalise<-phyloseq_to_metagenomeSeq(physeq_bac_filter_norm)
p_biom<-cumNormStatFast(physeq_bac_normalise)
biom_quant<-cumNorm(physeq_bac_normalise, p=p_biom)
biom_quant
normFactors(biom_quant)
physeq_bac_normalise<-MRcounts(biom_quant, norm=T)
head(physeq_bac_normalise)
physeq_bac_normalise
#create physeq object with normalized otu table
otu_table(physeq_bac_filter_norm) <- otu_table(physeq_bac_normalise,taxa_are_rows=T)
otu_table(physeq_bac_filter_norm)

#physeq_obj_ITS_uparse_R1_mSeq <- physeq_obj_ITS_uparse_R1_clean
#otu_table(physeq_obj_ITS_uparse_R1_mSeq) <- otu_table(biom_ITS_soil, taxa_are_rows=TRUE)

physeq_bac_filter_norm
head(otu_table(physeq_bac_filter_norm))
head(tax_table(physeq_bac_filter_norm))
head(sample_data(physeq_bac_filter_norm))

write.csv(otu_table(physeq_bac_filter_norm), file = "filtered2_otus.csv")




                 #######################################################################################
                 ##################           VISUALISATION             ################################
                 #######################################################################################
                       #https://joey711.github.io/phyloseq/plot_ordination-examples.html


                    ######################   UNCONSTRANED ORDINATION  ##########################
                     #######################################################################################
                 #https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#unconstrained_ordinations

#PCOA

#One of the best exploratory analyses for amplicon data is unconstrained ordinations. 
#Here we will look at ordinations of our full community samples. We will use the scale_reads() function in miseqR.R to scale to the smallest library size, which is the default. 
#If you want to scale to another depth, you can do so by setting the "n" argument

theme_set(theme_bw())
# Scale reads to even depth 

physeq_bac_scale <- transform_sample_counts(physeq_bac_filter_norm, function(x) 1E6 * x/sum(x) ) 
#Keep only the most abundant fifteen phyla.
#Genus.sum = tapply(taxa_sums(physeq_bac_scale), tax_table(physeq_bac_scale)[, "Genus"], sum, na.rm=TRUE)
#top19taxa = names(sort(Genus.sum, TRUE))[1:20]
#physeq_bac_scale = prune_taxa((tax_table(physeq_bac_scale)[, "Genus"] %in% top19taxa), physeq_bac_scale)

#physeq_bac_scale.ord <- ordinate(physeq_bac_scale, "NMDS", "bray")
#p1 = plot_ordination(physeq_bac_scale, physeq_bac_scale.ord, type="taxa", color="Genus", title="taxa")
#print(p1)

# Ordinate
ordination_pcoa <- ordinate(
  physeq = physeq_bac_scale, 
  method = "PCoA", 
  distance = "bray")#,weighted=TRUE)

# Plot 
colors <- c("Field2" = "#000000", "Field4" = "#FF0000", "Field5" = "#56B4E9", "Field6" = "#009E73", "Field8" = "#F0E442", "Field10" = "#0072B2", "Field13" = "#D55E00", "Field14" = "#CC79A7", "Field15" = "#999999",
            "Greenhouse2" = "#666666", "Greenhouse4" ="#A6761D", "Greenhouse5" = "#E6AB02", "Greenhouse6" = "#A6D854", "Greenhouse8" = "#E7298A", "Greenhouse10" = "#7570B3",
            "Greenhouse13" ="#A6CEE3", "Greenhouse14" = "#1F78B4", "Greenhouse15" = "#33A02C")


p <- plot_ordination(
  physeq = physeq_bac_scale,
  ordination = ordination_pcoa,
  color = "Population",
  shape = "Soil",  title = "PCoA of Bacteria Communities"
  ) 
p = p + geom_point(aes(colour = factor(Population)))
p = p + scale_colour_manual(values = colors)
p = p + geom_point(size=7, alpha=10)
#p = p + stat_ellipse(aes(group=Population), type="norm", alpha=0.7, linetype = 3, show.legend = FALSE)

print(p)                     


#scale_color_manual(values = c("#a65628", "red", "#ffae19","#4daf4a", "#1919ff", "darkorchid3", "magenta")

#p <- plot_ordination(physeq_bac_scale, ordu, color="Field", shape="Landscape")
#p <- p + geom_point(size=7, alpha=.7)
#p <- p + scale_colour_brewer(type="qual", palette="Set1")
#p <- p + ggtitle("MDS/PCoA on unweighted-UniFrac distance")

#print(p)


p <- ordinate(physeq_bac_scale, "PCoA", "unifrac", weighted=TRUE)
p <- plot_ordination(physeq_bac_scale, ordu, color="Field", shape="Landscape")
p <- p + geom_point(size=7, alpha=.7)
p <- p + scale_colour_brewer(type="qual", palette="Set1")
p <- p + ggtitle("MDS/PCoA on unweighted-UniFrac distance")

print(p)












#Indicator species
# Write out your phyloseq OTU table and export it
write.csv(biom_16S_uparse_filt@otu_table,'otus_ITS_out.csv')

#Import phyloseq OTU table as an OTU table/dataframe
SpOTU<-read.csv('otus_ITS_out.csv')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) # renames columns
SpOTUFlip<- SpOTUFlip[-1, ] #removes first row
SpOTUFlip_num<-as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$SampleID<-row.names(SpOTUFlip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata<-read.csv("16S_meta.csv") #read in metadata
head(metadata) # check

## Join based on SampleID
SpOTU_Final<-left_join(SpOTUFlip_num, metadata, by = c("SampleID" = "SampleID")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus = SpOTU_Final[,1:822] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat = SpOTU_Final$Population #the metadata column group you care about

SPind=multipatt(x=SPotus, cluster=SPwat, func = "r.g", control = how(nperm=9999))

summary(SPind, indvalcomp=TRUE)



                               ##############        STATISTICAL ANALYSIS         ###############################
                               ##################################################################################
                                              #https://microbiome.github.io/tutorials/#



                                         ########### PERMANOVA metagenomeseq ################
                             ###################################################################################
library("vegan")
library("RVAideMemoire")

options(scipen = 999)
set.seed(1)

#Subset agricultural fields for analysis

physeq_bac_filter_norm_agric <- subset_samples(physeq_bac_filter_norm, LAND_USE=="Agricultural")


physeq_bac_filter_norm_metadata_agric <- as.data.frame(as.matrix(sample_data(physeq_bac_filter_norm_agric)))
head(physeq_bac_filter_norm_metadata_agric)

model.matrix(~ ST * Mhapla, data = physeq_bac_filter_norm_metadata_agric)
model.matrix(~ ST + Mhapla + ST:Mhapla, physeq_bac_filter_norm_metadata_agric)

# Adonis test
adonis(t(otu_table(physeq_bac_filter_norm_agric)) ~ ST * Mhapla, data=physeq_bac_filter_norm_metadata_agric, permutations=9999) # by = "margin"
adonis(t(otu_table(physeq_bac_filter_norm_agric)) ~ ST + Mhapla + ST:Mhapla, data=physeq_bac_filter_norm_metadata_agric, permutations=9999) 
adonis(t(otu_table(physeq_bac_filter_norm_agric)) ~ ST + Mhapla + ST:Mhapla, data=physeq_bac_filter_norm_metadata_agric, permutations=9999)

write.csv(otu_table(physeq_bac_filter_norm_agric), file = "otu_check.csv")
obj1_otu <- as.data.frame(t(otu_table(physeq_bac_filter_norm_agric)))
head(obj1_otu)          

# Homogeneity of dispersion test
vegan::vegdist(obj1_otu, method="bray") -> dist_otu_physeq_bac_filter_norm_agric

permdisp_otu_ST<- betadisper(dist_otu_physeq_bac_filter_norm_agric, physeq_bac_filter_norm_metadata_agric$ST, type = "centroid")
permdisp_otu_Mhapla <- betadisper(dist_otu_physeq_bac_filter_norm_agric, physeq_bac_filter_norm_metadata_agric$Mhapla, type = "centroid")
permdisp_otu_ST
permdisp_otu_Mhapla

anova(permdisp_otu_ST, permutations=9999)
anova(permdisp_otu_Mhapla, permutations=9999)

permutest(permdisp_otu_prok_soil_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_soil_M)
plot(TukeyHSD(permdisp_otu_prok_soil_M), las=1)
boxplot(permdisp_otu_prok_soil_M)























# alpha diversity plots
### alpha diversity
# making graph following reviewers reccomendations
library("gridExtra")
library("grid")
library("cowplot")

alpha_supp <- biom_16S_uparse_filt
otu_prok <- as.data.frame(otu_table(alpha_supp ))

otu_prok
meta_prok <- as.data.frame(sample_data(alpha_supp ))
alpha_prok <- meta_prok
alpha_prok
alpha_prok$readNO <- sample_sums(alpha_supp)
alpha_prok$Observed <- specnumber(otu_prok, MARGIN = 2)
alpha_prok$Shannon <- diversity(otu_prok, index="shannon", MARGIN = 2)
jevenness_prok <- diversityresult(t(otu_prok), method = "each site", index = "Jevenness")
alpha_prok$Jevenness <- jevenness_prok$Jevenness
#alpha_prok <- alpha_prok[order(alpha_prok$ReadNO), ]
alpha_prok

alpha_prok  
alpha_prok$alpha_label <- factor(alpha_prok$LAND_USE,
                                 level=c("Agricultural","Forested"))
p <- ggplot(alpha_prok, aes(x=alpha_label, y=readNO)) + 
  theme_classic()+
  #scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      #values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  geom_point(size = 2, shape = 16) +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold"))+
  geom_boxplot()
p

label_names <- c(Observed="Richness", Shannon="Shannon")
label_names


ps.noncontam_alpha_soil <- biom_16S_uparse_filt
sample_data(ps.noncontam_alpha_soil)$alpha_label <- factor(sample_data(ps.noncontam_alpha_soil)$alpha_label,
                                                           level=c("Agricultural","Forested"))
estimate_richness(ps.noncontam_alpha_soil, split = TRUE, measures = NULL)
alpha_soil_prok = plot_richness(ps.noncontam_alpha_soil, x= "alpha_label", 
                                , measures = c( "Observed")) +
  
  ylim(0,9000) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Observed OTUs") +
  scale_x_discrete("Sample", labels = c("Agricultural" = "Forested")) +
  #scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      #values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  
  theme(axis.text.x = element_text(angle = 90,vjust =1.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  
  theme_set(theme_classic())+
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  
  theme(legend.title=element_blank())
plot(alpha_soil_prok)

alpha_soil_prok_shan = plot_richness(ps.noncontam_alpha_soil, x= "alpha_label", 
                                      measures = c( "Shannon")) +
  
  ylim(1,8) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  
  theme(legend.title=element_blank())
plot(alpha_soil_prok_shan)



























# create vegan objects 
otu_prokaryote_out <- as.data.frame(otu_table(biom_16S_uparse_filt))
taxa_prokaryote_out <- as.data.frame(as.matrix(tax_table(biom_16S_uparse_filt)))
metadata_prokaryote_out <- as.data.frame(as.matrix(sample_data(biom_16S_uparse_filt)))
dim(otu_prokaryote_out)






























# >>> RAREFACTION CURVES -----------------------------------------------------------------------------

# *** FIGURE 1 rarecurve prokaryotes  -----------------------------------------------------------
rarecurve(t(otu_prokaryote_out), col = metadata_prokaryote_out$Field, label = FALSE, sample=min(colSums(otu_prokaryote_out)), step = 50,
          main="Bacteria", ylab = "Number of OTUs", xlab = "Number of DNA reads", cex=0.6) -> rare_prokaryote
legend("bottomright", legend=c("R", "S", "SCSC","RCRC"),
       col=c("red", "black",  "red","blue"), lty=1, cex=0.8, box.lty=1) #box.lty=0 remove the legend border


#>>>> INDICATOR SPECIES FROM SOIL MICROBIOME

#Creating a phyloseq with only 25 otus
#read 25 otus from .txt file
library(utils)
library(base)
??read.table

#Subset 9 classes associated with M. hapla presence in field soils NB at the class level only seven classes are present
rownames_isa_fdr <- read.table('indicator_soil_otus_with_their_taxonomy.txt', header = TRUE)

my_subset = subset_taxa(biom_16S_uparse, Class=="Ellin6529" | Class=="MB-A2-108" | Class=="Actinobacteria"| Class=="Thermoleophilia"| Class=="Planctomycetia"| Class=="Gitt-GS-136"| Class=="Alphaproteobacteria" )



#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(my_subset, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(my_subset)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)
#write.table(physeq_bacRm,"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_Class <- physeq_bacRm  %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Class)           # Sort data frame alphabetically by Class
pm_Class

dat_Soil_Pre_Abs_bp <- data.table(pm_Class)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Class:= "Other"]

#manually ordering levels
dat_Soil_Pre_Abs_bp1 <- dat_Soil_Pre_Abs_bp
dat_Soil_Pre_Abs_bp1
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13","Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")

#pm_Class <- get_top_taxa(dat_Soil_Pre_Abs_bp, 20, relative = TRUE, other_label = "Other")
#Plot overall abundance in ggplot 
p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Class))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C",
                               "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" ,"#FBB4AE", "#B3CDE3" ,"#CCEBC5", "#DECBE4", "#FED9A6",
                               "#FFFFCC" ,"#E5D8BD", "#FDDAEC", "#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
                               "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999" ,"#66C2A5",
                               "#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072",
                               "#80B1D3", "#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F")) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=21)) + scale_x_discrete(limits = positions)







#####################Only actinomycetes
rownames_isa_fdr <- read.table('indicator_soil_otus_with_their_taxonomy.txt', header = TRUE)

my_subset2 = subset_taxa(biom_16S_uparse, Class=="Actinobacteria" )


#Relative abundance fields otu table to 100%
physeq_bacRm = merge_samples(my_subset2, "Population")
sample_data(physeq_bacRm)$ Population <- levels(sample_data(my_subset2)$Population)
physeq_bacRm = transform_sample_counts(physeq_bacRm, function(x) 100 * x/sum(x))
View(physeq_bacRm)
#write.table(physeq_bacRm,"physeq_bacRm_relabu.csv") #output the rel. abu table for further analysis in excel
#Plot relative abundance
pm_Class <- physeq_bacRm  %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate taxa at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()   %>%                                            # Melt phyloseq object to long format for producing graphics with ggplot2
  arrange(Order)           # Sort data frame alphabetically by Order
pm_Class

dat_Soil_Pre_Abs_bp <- data.table(pm_Class)
dat_Soil_Pre_Abs_bp [(Abundance <= 0.00), Order:= "Other"]

#manually ordering levels
dat_Soil_Pre_Abs_bp1 <- dat_Soil_Pre_Abs_bp
dat_Soil_Pre_Abs_bp1
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13","Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
positions <- c("Greenhouse4","Greenhouse5","Greenhouse6","Greenhouse10","Greenhouse14","Greenhouse15","Greenhouse2","Greenhouse8","Greenhouse13")
#positions <- c("Field4","Field5","Field6","Field10","Field14","Field15","Field2","Field8","Field13")

#pm_Class <- get_top_taxa(dat_Soil_Pre_Abs_bp, 20, relative = TRUE, other_label = "Other")
#Plot overall abundance in ggplot 
p <- ggplot(data=dat_Soil_Pre_Abs_bp, aes(x=Sample, y=Abundance, fill=Order))
p + geom_bar(aes(), stat="identity", position="stack") + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  geom_bar(stat = "identity", width = 0.1) +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086" ,"#386CB0" ,"#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3",
                               "#E7298A", "#66A61E", "#E6AB02", "#A6761D" ,"#666666" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C",
                               "#FF7F00" ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928" ,"#FBB4AE", "#B3CDE3" ,"#CCEBC5", "#DECBE4", "#FED9A6",
                               "#FFFFCC" ,"#E5D8BD", "#FDDAEC", "#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8", "#F4CAE4" ,"#E6F5C9" ,"#FFF2AE", "#F1E2CC",
                               "#CCCCCC", "#E41A1C" ,"#377EB8" ,"#4DAF4A" ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF" ,"#999999" ,"#66C2A5",
                               "#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#A6D854" ,"#FFD92F" ,"#E5C494" ,"#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA" ,"#FB8072",
                               "#80B1D3", "#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD" ,"#CCEBC5", "#FFED6F")) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=21)) + scale_x_discrete(limits = positions)












#>>> BETA DIVERSITY ------------------------------------------------------------------------------
library("ggrepel")
library("ggplot2")
colnames(sample_data(physeq_fungi_mSeq))[2] <- "Matrix" #set the column names of matrix
sample_data(physeq_fungi_mSeq) <- sample_data(physeq_fungi_mSeq)[,c(1,2,3,4,5,6,7,8)]
sample_data(physeq_fungi_mSeq)

pcoa_fungi_out = phyloseq::ordinate(physeq_fungi_mSeq, method ="PCoA", distance="bray")

p_pcoa_ITS_out = plot_ordination(physeq_fungi_mSeq, pcoa_fungi_out, color="Rotation", shape ="SCN") + 
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=3, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values=c("darkgoldenrod1", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon",
                               "dodgerblue3", "steelblue1", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkgrey","red", "grey", "seagreen3")) +
  scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15)) +
  stat_ellipse(aes(group=Rotation), type="norm", alpha=0.8, linetype = 3, show.legend = FALSE) +
  #geom_text_repel(aes(label=sample_data(physeq_fungi_mSeq)$Description), size = 2) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) +
  #theme(legend.title=element_blank())
  theme(legend.position="bottom") 
p_pcoa_ITS_out
p_pcoa_ITS_out + stat_ellipse(geom = "polygon", level=0.70, type="norm", alpha=0.04)


colnames(sample_data(physeq_prokaryote_mSeq))[2] <- "Matrix" #set the column names of matrix
sample_data(physeq_prokaryote_mSeq) <- sample_data(physeq_prokaryote_mSeq)[,c(1,2,3,4,5,6,7,8)]
sample_data(physeq_prokaryote_mSeq)

pcoa_prokaryote_out = phyloseq::ordinate(physeq_prokaryote_mSeq, method ="PCoA", distance="bray")

p_pcoa_16S_out = plot_ordination(physeq_prokaryote_mSeq, pcoa_prokaryote_out, color="Rotation", shape ="SCN") + 
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=3, alpha=1) + # ,aes(shape=Age))
  scale_colour_manual(values=c("darkgoldenrod1", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon",
                               "dodgerblue3", "steelblue1", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkgrey","red", "grey", "seagreen3")) +
  scale_shape_manual(values=c(1, 2,3,4,5,6,7,8,9,10,11,12,13,14,15)) +
  #stat_ellipse(aes(group=Rotation), type="norm", alpha=0.8, linetype = 2, show.legend = FALSE) +
  #geom_text_repel(aes(label=sample_data(physeq_prokaryote_mSeq)$Description), size = 2) + 
  theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 8)) +
  #theme(legend.title=element_blank())
  theme(legend.position="bottom") 
p_pcoa_16S_out




# *** FIGURE 3 - ordinations ---------------------------------------------------------------------
library("ggpubr")

ggarrange(p_pcoa_16S_out,
          p_pcoa_ITS_out,
          labels = c("A", "B"),
          widths = c(1,1),
          align = "none", 
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = c("right"))


# >>> PERMANOVA ----------------------------------------------------------------------------------
library (vegan)
# fungal communities 
vegan::vegdist(t(otu_fungi_out_mSeq), method="bray") -> dist_fungi

model.matrix(~ Tillage*SCN*Rotation, data=metadata_fungi_out_mSeq)
adonis(dist_fungi ~ Tillage*SCN*Rotation, data=metadata_fungi_out_mSeq, permutations=999) -> adonis_fungi
adonis_fungi

densityplot(permustats(adonis_fungi), main=list(label="Fungi var. Tillage/SCN/Rotation", cex=1)) -> density_fungi
density_fungi

betadisper(dist_fungi, metadata_fungi_out_mSeq$Tillage) -> betadisper_fungi_tillage
anova(betadisper_fungi_tillage, permutations = 9999)
permutest(betadisper_fungi_tillage, permutations = 9999, pairwise = T)
plot(betadisper_fungi_tillage)
plot(TukeyHSD(betadisper_fungi_tillage), las=0)
boxplot(betadisper_fungi_tillage)

betadisper(dist_fungi, metadata_fungi_out_mSeq$SCN) -> betadisper_fungi_SCN
anova(betadisper_fungi_SCN, permutations = 9999)
permutest(betadisper_fungi_SCN, permutations = 9999, pairwise = T)
plot(betadisper_fungi_SCN)
plot(TukeyHSD(betadisper_fungi_SCN), las=0)
boxplot(betadisper_fungi_SCN)

betadisper(dist_fungi, metadata_fungi_out_mSeq$Rotation) -> betadisper_fungi_Rotation
anova(betadisper_fungi_Rotation, permutations = 9999)
permutest(betadisper_fungi_Rotation, permutations = 9999, pairwise = T)
plot(betadisper_fungi_Rotation)
#plot(betadisper_fungi_Rotation,hull = FALSE, ellipse = TRUE)
plot(TukeyHSD(betadisper_fungi_Rotation), las=0)
boxplot(betadisper_fungi_Rotation)

# prokaryotic communities 
vegan::vegdist(t(otu_prokaryote_out_mSeq), method="bray") -> dist_prokaryote

model.matrix(~ Tillage*SCN*Rotation, data=metadata_prokaryote_out_mSeq)
adonis(dist_prokaryote ~ Tillage*SCN*Rotation, data=metadata_prokaryote_out_mSeq, permutations=999) -> adonis_prokaryote
adonis_prokaryote
densityplot(permustats(adonis_prokaryote), main=list(label="Prokaryote var. Tillage/SCN/Rotation", cex=1)) -> density_prokaryote
density_prokaryote

betadisper(dist_prokaryote, metadata_prokaryote_out_mSeq$Tillage) -> betadisper_prokaryote_Tillage
anova(betadisper_prokaryote_Tillage, permutations = 9999)
permutest(betadisper_prokaryote_Tillage, permutations = 9999, pairwise = T)
plot(betadisper_prokaryote_Tillage)
plot(TukeyHSD(betadisper_prokaryote_Tillage), las=0)
boxplot(betadisper_prokaryote_Tillage)

betadisper(dist_prokaryote, metadata_prokaryote_out_mSeq$SCN) -> betadisper_prokaryote_SCN
anova(betadisper_prokaryote_SCN, permutations = 9999)
permutest(betadisper_prokaryote_SCN, permutations = 9999, pairwise = T)
plot(betadisper_prokaryote_SCN)
plot(TukeyHSD(betadisper_prokaryote_SCN), las=0)
boxplot(betadisper_prokaryote_SCN)

betadisper(dist_prokaryote, metadata_prokaryote_out_mSeq$Rotation) -> betadisper_prokaryote_Rotation
anova(betadisper_prokaryote_Rotation, permutations = 9999)
permutest(betadisper_prokaryote_Rotation, permutations = 9999, pairwise = T)
plot(betadisper_prokaryote_Rotation)
plot(TukeyHSD(betadisper_prokaryote_Rotation), las=0)
boxplot(betadisper_prokaryote_Rotation)


# Figure S4 - betadisper -------------------------------------------------------------------------
par(mfrow=c(3,3))
plot(betadisper_prokaryote_Tillage, main="(A) Prokaryotes\n\nPCoA (Tillage)", las=1) #, cex.main=1
boxplot(betadisper_prokaryote_Tillage, main="\nBoxplot (Tillage)", xlab="Tillage")
plot(TukeyHSD(betadisper_prokaryote_Tillage), las=0)

plot(betadisper_prokaryote_SCN, main="(B)Prokaryotes\nPCoA (SCN)", las=1) #, cex.main=1
boxplot(betadisper_prokaryote_SCN, main="\nboxplot (SCN)", xlab="Stage")
plot(TukeyHSD(betadisper_prokaryote_SCN), las=0)

plot(betadisper_prokaryote_Rotation, main="(C)Prokaryotes\nPCoA (Rotation)", las=1) #, cex.main=1
boxplot(betadisper_prokaryote_Rotation, main="\nboxplot (Rotation)", xlab="Rotation")
plot(TukeyHSD(betadisper_prokaryote_Rotation), las=0)



plot(betadisper_fungi_tillage, main="(A) Fungi\n\nPCoA (Tillage)", las=1) #, cex.main=1
boxplot(betadisper_fungi_tillage, main="\nBoxplot (Tillage)", xlab="Tillage")
plot(TukeyHSD(betadisper_fungi_tillage), las=0)

plot(betadisper_fungi_SCN, main="(B)Fungi\nPCoA (SCN)", las=1) #, cex.main=1
boxplot(betadisper_fungi_SCN, main="\nboxplot (SCN)", xlab="Stage")
plot(TukeyHSD(betadisper_fungi_SCN), las=0)

plot(betadisper_fungi_Rotation, main="(C)Fungi\nPCoA (Rotation)", las=1) #, cex.main=1
boxplot(betadisper_fungi_Rotation, main="\nboxplot (Rotation)", xlab="Rotation")
plot(TukeyHSD(betadisper_fungi_Rotation), las=0)

dev.off()


# >>> ALPHA DIVERSITY ----------------------------------------------------------------------------
library("BiodiversityR")
library(vegan)
dim(otu_fungi_out)
dim(metadata_fungi_out)
count(metadata_fungi_out, vars = Tillage)
count(metadata_fungi_out, vars = SCN)
count(metadata_fungi_out, vars = Rotation)

alpha_div_fungi <- metadata_fungi_out[,c(5:7)]
alpha_div_fungi$readNO <- sample_sums(physeq_fungi)
alpha_div_fungi$Observed <- specnumber(otu_fungi_out, MARGIN = 2)
alpha_div_fungi$Rarefied <- rarefy(otu_fungi_out,sample=min(alpha_div_fungi$readNO), MARGIN = 2)
alpha_div_fungi$Shannon <- diversity(otu_fungi_out, index="shannon")
#jevenness_fungi <- diversityresult(t(otu_fungi_out), method = "each site", index = "Jevenness")
#alpha_div_fungi$Jevenness <- jevenness_fungi$Jevenness
alpha_div_fungi <- alpha_div_fungi[order(alpha_div_fungi$readNO), ]
alpha_div_fungi

??diversityresult

# get descriptive stats
library("psych")
describeBy(alpha_div_fungi, alpha_div_fungi$Tillage)
describeBy(alpha_div_fungi, alpha_div_fungi$SCN)
describeBy(alpha_div_fungi, alpha_div_fungi$Rotation)

# check homogenity of variances 
# check homogenity of variances 
fligner.test(Observed ~ Rotation, data=alpha_div_fungi)
fligner.test(shannon ~ Rotation, data=alpha_div_fungi)
fligner.test(Jevenness ~ Rotation, data=alpha_div_fungi)

# get significant differences
aov_fungi_rich <- aov(Observed ~ Rotation, data=alpha_div_fungi)
summary(aov_fungi_rich)
aov_fungi_shan <- aov(shannon ~ Rotation, data=alpha_div_fungi)
summary(aov_fungi_shan)
aov_fungi_Jeven <- aov(Jevenness ~ Rotation, data=alpha_div_fungi)
summary(aov_fungi_Jeven)

??count
library(phyloseq)
library(ggplot2)## for plotting
library(magrittr)
library(ggpubr)##for combining the plots
library(vegan)##for community ecology based codes
library(limma)
library(edgeR)
library(Hmisc)
library(igraph)
library(labdsv)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(microbiome)
library(knitr)


#Prokaryote
library("BiodiversityR")
dim(otu_prokaryote_out)
dim(metadata_prokaryote_out)
count(metadata_prokaryote_out, vars = Agricultural)
count(metadata_prokaryote_out, vars = SCN)
count(metadata_prokaryote_out, vars = Rotation)
?count

alpha_div_prokaryote <- metadata_prokaryote_out[,c(5:7)]
alpha_div_prokaryote$readNO <- sample_sums(physeq_prokaryote)
alpha_div_fungi$Observed <- specnumber(otu_fungi_out, MARGIN = 2)
alpha_div_fungi$Rarefied <- rarefy(otu_fungi_out,sample=min(alpha_div_fungi$readNO), MARGIN = 2)
alpha_div_fungi$Shannon <- diversity(otu_fungi_out, index="shannon")
jevenness_fungi <- diversityresult(t(otu_fungi_out), method = "each site", index = "Jevenness")
alpha_div_fungi$Jevenness <- jevenness_fungi$Jevenness
alpha_div_fungi <- alpha_div_fungi[order(sums_fungi$ReadNO), ]
alpha_div_fungi


# get descriptive stats
library("psych")
describeBy(alpha_div_prokaryote, alpha_div_prokaryote$Tillage)
describeBy(alpha_div_prokaryote, alpha_div_prokaryote$SCN)
describeBy(alpha_div_prokaryote, alpha_div_prokaryote$Rotation)

# check homogenity of variances 
# check homogenity of variances 
fligner.test(Observed ~ Rotation, data=alpha_div_prokaryote)
fligner.test(Shannon ~ Rotation, data=alpha_div_prokaryote)
fligner.test(Jevenness ~ Rotation, data=alpha_div_prokaryote)

# get significant differences
aov_prokaryote_rich <- aov(Observed ~ Rotation, data=alpha_div_prokaryote)
summary(aov_prokaryote_rich)
aov_prokaryote_shan <- aov(Shannon ~ Rotation, data=alpha_div_prokaryote)
summary(aov_prokaryote_shan)
aov_prokaryote_Jeven <- aov(Jevenness ~ Rotation, data=alpha_div_prokaryote)
summary(aov_prokaryote_Jeven)

#dim(otu_prokaryote_out)
#dim(metadata_prokaryote_out)
#count(metadata_prokaryote_out, vars = Tillage)
#count(metadata_prokaryote_out, vars = SCN)
#count(metadata_prokaryote_out, vars = Rotation)

#alpha_div_prokaryote <- metadata_prokaryote_out[,c(1,4,6:7)]
#alpha_div_prokaryote$readNO <- sample_sums(physeq_prokaryote)
#alpha_div_prokaryote$Observed <- specnumber(otu_prokaryote_out, MARGIN = 2)
#alpha_div_prokaryote$Rarefied <- rarefy(otu_prokaryote_out,sample=min(alpha_div_prokaryote$readNO), MARGIN = 2)
#alpha_div_prokaryote$Shannon <- diversity(otu_prokaryote_out, index="shannon", MARGIN = 2)
#jevenness_prokaryote <- diversityresult(t(otu_prokaryote_out), method = "each site", index = "Jevenness")
#alpha_div_prokaryote$Jevenness <- jevenness_prokaryote$Jevenness
#alpha_div_prokaryote <- alpha_div_prokaryote[order(sums_prokaryote$ReadNO), ]
#alpha_div_prokaryote

#describeBy(alpha_div_prokaryote, alpha_div_prokaryote$Stage)
#describeBy(alpha_div_prokaryote, alpha_div_prokaryote$Origin)

# check homogenity of variances 
#fligner.test(Observed ~ Origin, data=alpha_div_prokaryote)
#fligner.test(Shannon ~ Origin, data=alpha_div_prokaryote)
#fligner.test(Jevenness ~ Origin, data=alpha_div_prokaryote)

# get significant differences
#library("agricolae")

#aov_prokaryote_rich <- aov(Observed ~ Stage, data=alpha_div_prokaryote)
#summary(aov_prokaryote_rich)
#aov_prokaryote_shan <- aov(Shannon ~ Stage, data=alpha_div_prokaryote)
#summary(aov_prokaryote_shan)
#aov_prokaryote_Jeven <- aov(Jevenness ~ Stage, data=alpha_div_prokaryote)
#summary(aov_prokaryote_Jeven)

#aov_prokaryote_rich <- aov(Observed ~ Origin, data=alpha_div_prokaryote)
#summary(aov_prokaryote_rich)
#HSD.test(aov_prokaryote_rich, "Origin") -> tukeyHSD_prokaryote_rich
#tukeyHSD_prokaryote_rich

#aov_prokaryote_shan <- aov(Shannon ~ Origin, data=alpha_div_prokaryote)
#summary(aov_prokaryote_shan)
#HSD.test(aov_prokaryote_shan, "Origin") -> tukeyHSD_prokaryote_shan
#tukeyHSD_prokaryote_shan

#aov_prokaryote_Jeven <- aov(Jevenness ~ Origin, data=alpha_div_prokaryote)
#summary(aov_prokaryote_Jeven)
#HSD.test(aov_prokaryote_Jeven, "Origin") -> tukeyHSD_prokaryote_Jeven
#tukeyHSD_prokaryote_Jeven


# >>> INDICATOR SPECIES ANALYSIS (ISA) -----------------------------------------------------------
library("indicspecies")
library("phyloseq")

#Preparing data from phyloseq for indispecies
#Subset Degraded_Min_no_Mh 
#isa_fungi_1<- subset_samples(biom_16S_uparse_filt, Mhapla=="Positive")

#isa_fungi_1 <- subset_samples(physeq_bac_filter, Mhapla=="Positive")
#Disturbed_Min_with_Mh
#Disturbed_Min_no_Mh
#Disturbed_Muc_with_Mh
#Degraded_Muc_with_Mh
#Maturing_Min_no_Mh


# Write out your phyloseq OTU table and export it
write.csv(biom_16S_uparse_filt@otu_table,'otus_ITS_out.csv')

#Import phyloseq OTU table as an OTU table/dataframe
SpOTU<-read.csv('otus_ITS_out.csv')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) # renames columns
SpOTUFlip<- SpOTUFlip[-1, ] #removes first row
SpOTUFlip_num<-as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$SampleID<-row.names(SpOTUFlip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata<-read.csv("ITS_metadata.csv") #read in metadata
head(metadata) # check

## Join based on SampleID
SpOTU_Final<-left_join(SpOTUFlip_num, metadata, by = c("SampleID" = "SampleID")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus = SpOTU_Final[,1:1377] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat = SpOTU_Final$Soilfoodweb #the metadata column group you care about

SPind=multipatt(x=SPotus, cluster=SPwat, func = "r.g", control = how(nperm=9999))

summary(SPind, indvalcomp=TRUE)





#disturbed muck with Mhapla
isa_fungi_1 <- multipatt(as.data.frame(t(physeq_bac_filter$otutable), ITS_metadata_fungi$Indicator_group, control=how(nperm=9999))
summary(isa_fungi_T, indvalcomp=TRUE)


isa_fungi_T -> isa_fungi_fdr_T
isa_fungi_fdr_T$sign$p.value<-p.adjust(isa_fungi_T$sign$p.value, "fdr")
summary(isa_fungi_fdr_T)

as.data.frame(t(physeq_bac_filter$otu_table)



# >>> INDICATOR SPECIES ANALYSIS (ISA) -----------------------------------------------------------
library("indicspecies")
isa_fungi_T <- multipatt(as.data.frame(t(otu_fungi_out)), metadata_fungi_out$Tillage, control=how(nperm=9999))
summary(isa_fungi_T, indvalcomp=TRUE)
isa_fungi_T -> isa_fungi_fdr_T
isa_fungi_fdr_T$sign$p.value<-p.adjust(isa_fungi_T$sign$p.value, "fdr")
summary(isa_fungi_fdr_T)

isa_fungi_S <- multipatt(as.data.frame(t(otu_fungi_out)), metadata_fungi_out$SCN, control=how(nperm=9999))
summary(isa_fungi_S, indvalcomp=TRUE)
isa_fungi_S -> isa_fungi_fdr_S
isa_fungi_fdr_S$sign$p.value<-p.adjust(isa_fungi_S$sign$p.value, "fdr")
summary(isa_fungi_fdr_S)

isa_fungi_R <- multipatt(as.data.frame(t(otu_fungi_out)), metadata_fungi_out$Rotation, control=how(nperm=9999))
summary(isa_fungi_R, indvalcomp=TRUE)
isa_fungi_R -> isa_fungi_fdr_R
isa_fungi_fdr_R$sign$p.value<-p.adjust(isa_fungi_R$sign$p.value, "fdr")
summary(isa_fungi_fdr_R)


isa_prokaryote_T <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$Tillage, control=how(nperm=9999))
summary(isa_prokaryote_T, indvalcomp=TRUE)
isa_prokaryote_T -> isa_prokaryote_fdr_T
isa_prokaryote_fdr_T$sign$p.value<-p.adjust(isa_prokaryote_T$sign$p.value, "fdr")
summary(isa_prokaryote_fdr_T)

isa_prokaryote_S <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$SCN, control=how(nperm=9999))
summary(isa_prokaryote_S, indvalcomp=TRUE)
isa_prokaryote_S -> isa_prokaryote_fdr_S
isa_prokaryote_fdr_S$sign$p.value<-p.adjust(isa_prokaryote_S$sign$p.value, "fdr")
summary(isa_prokaryote_fdr_S)

isa_prokaryote_R <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$Rotation, control=how(nperm=9999))
summary(isa_prokaryote_R, indvalcomp=TRUE)
isa_prokaryote_R -> isa_prokaryote_fdr_R
isa_prokaryote_fdr_R$sign$p.value<-p.adjust(isa_prokaryote_R$sign$p.value, "fdr")
summary(isa_prokaryote_fdr_R)

#isa_prokaryote_ST <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$Stage, control=how(nperm=9999))
#summary(isa_prokaryote_ST, indvalcomp=TRUE)
#isa_prokaryote_ST -> isa_prokaryote_fdr_ST
#isa_prokaryote_fdr_ST$sign$p.value<-p.adjust(isa_prokaryote_ST$sign$p.value, "fdr")
#summary(isa_prokaryote_fdr_ST)

#isa_prokaryote_OR <- multipatt(as.data.frame(t(otu_prokaryote_out)), metadata_prokaryote_out$Origin, control=how(nperm=9999))
#summary(isa_prokaryote_OR, indvalcomp=TRUE)
#isa_prokaryote_OR -> isa_prokaryote_fdr_OR
#isa_prokaryote_fdr_OR$sign$p.value<-p.adjust(isa_prokaryote_OR$sign$p.value, "fdr")
#summary(isa_prokaryote_fdr_OR)

sink(file="isa_fungi_fdr_T.csv") 
summary(isa_fungi_fdr_T)
sink(file="isa_fungi_fdr_S.csv") 
summary(isa_fungi_fdr_S)
sink(file="isa_fungi_fdr_R.csv") 
summary(isa_fungi_fdr_R)

sink(file="isa_prokaryote_fdr_T.csv") 
summary(isa_prokaryote_fdr_T)
sink(file="isa_prokaryote_fdr_S.csv") 
summary(isa_prokaryote_fdr_S)
sink(file="isa_prokaryote_fdr_R.csv") 
summary(isa_prokaryote_fdr_R)

#sink()

# extracting ISA OTUs ----------------------------------------------------------------------------
#####Only TILLAGE
result_isa_fungi_fdr_T <- isa_fungi_fdr_T$sign[which(isa_fungi_fdr_T$sign$p.value <= 0.05), ]
head(result_isa_fungi_fdr_T)
dim(result_isa_fungi_fdr_T)

result_isa_fungi_fdr_T[result_isa_fungi_fdr_T$s.NT==0 &
                         result_isa_fungi_fdr_T$s.T==1 ,] -> T

result_isa_fungi_fdr_T[result_isa_fungi_fdr_T$s.NT==1 &
                         result_isa_fungi_fdr_T$s.T==0 ,] -> NT
isa_df_T <- rbind(NT, T)
dim(isa_df_T)
isa_df_T

#Only SCN
result_isa_fungi_fdr_S <- isa_fungi_fdr_S$sign[which(isa_fungi_fdr_S$sign$p.value <= 0.05), ]
head(result_isa_fungi_fdr_S)
dim(result_isa_fungi_fdr_S)

result_isa_fungi_fdr_T[result_isa_fungi_fdr_S$s.Yes==0 &
                         result_isa_fungi_fdr_S$s.No==1 ,] -> No_SCN

result_isa_fungi_fdr_T[result_isa_fungi_fdr_S$s.Yes==1 &
                         result_isa_fungi_fdr_S$s.No==0 ,] -> Yes_SCN

isa_df_S <- rbind(Yes_SCN, No_SCN)
#isa_df <- rbind(NT, T, Yes_SCN, No_SCN, C, R, S, RCRC, SCSC, C_SCSC, S_SCSC, C_RCRC_SCSC, R_S_SCSC,R_RCRC_S_SCSC)
dim(isa_df_S)
isa_df_S

#Only ROTATION
result_isa_fungi_fdr_R <- isa_fungi_fdr_R$sign[which(isa_fungi_fdr_R$sign$p.value <= 0.05), ]
head(result_isa_fungi_fdr_R)
dim(result_isa_fungi_fdr_R)

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==1 &
                               result_isa_fungi_fdr_R$s.S==0 &
                               result_isa_fungi_fdr_R$s.C==0 &
                               result_isa_fungi_fdr_R$s.SCSC==0 &
                               result_isa_fungi_fdr_R$s.RCRC==0 ,] -> R

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==0 &
                         result_isa_fungi_fdr_R$s.C==1 &
                         result_isa_fungi_fdr_R$s.SCSC==0 &
                         result_isa_fungi_fdr_R$s.RCRC==0 ,] -> C

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==1 &
                         result_isa_fungi_fdr_R$s.C==0 &
                         result_isa_fungi_fdr_R$s.SCSC==0 &
                         result_isa_fungi_fdr_R$s.RCRC==0 ,] -> S

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==0 &
                         result_isa_fungi_fdr_R$s.C==0 &
                         result_isa_fungi_fdr_R$s.SCSC==1 &
                         result_isa_fungi_fdr_R$s.RCRC==0 ,] -> SCSC

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==0 &
                         result_isa_fungi_fdr_R$s.C==0 &
                         result_isa_fungi_fdr_R$s.SCSC==0 &
                         result_isa_fungi_fdr_R$s.RCRC==1 ,] -> RCRC

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==0 &
                         result_isa_fungi_fdr_R$s.C==1 &
                         result_isa_fungi_fdr_R$s.SCSC==1 &
                         result_isa_fungi_fdr_R$s.RCRC==0 ,] -> C_SCSC

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==1 &
                         result_isa_fungi_fdr_R$s.C==0 &
                         result_isa_fungi_fdr_R$s.SCSC==1 &
                         result_isa_fungi_fdr_R$s.RCRC==0 ,] -> S_SCSC

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==0 &
                         result_isa_fungi_fdr_R$s.S==0 &
                         result_isa_fungi_fdr_R$s.C==1 &
                         result_isa_fungi_fdr_R$s.SCSC==1 &
                         result_isa_fungi_fdr_R$s.RCRC==1 ,] -> C_RCRC_SCSC

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==1 &
                         result_isa_fungi_fdr_R$s.S==1 &
                         result_isa_fungi_fdr_R$s.C==0 &
                         result_isa_fungi_fdr_R$s.SCSC==1 &
                         result_isa_fungi_fdr_R$s.RCRC==0 ,] -> R_S_SCSC

result_isa_fungi_fdr_R[result_isa_fungi_fdr_R$s.R==1 &
                         result_isa_fungi_fdr_R$s.S==0 &
                         result_isa_fungi_fdr_R$s.C==0 &
                         result_isa_fungi_fdr_R$s.SCSC==1 &
                         result_isa_fungi_fdr_R$s.RCRC==1 ,] -> R_RCRC_S_SCSC

isa_df_R <- rbind(R, C, S, SCSC, RCRC, C_SCSC, S_SCSC, C_RCRC_SCSC, R_S_SCSC,R_RCRC_S_SCSC)
#isa_df <- rbind(NT, T, Yes_SCN, No_SCN, C, R, S, RCRC, SCSC, C_SCSC, S_SCSC, C_RCRC_SCSC, R_S_SCSC,R_RCRC_S_SCSC)
dim(isa_df_R)
isa_df_R


######Only TILLAGE
result_isa_prokaryote_fdr_T <- isa_prokaryote_fdr_T$sign[which(isa_prokaryote_fdr_T$sign$p.value <= 0.05), ]
head(result_isa_prokaryote_fdr_T)
dim(result_isa_prokaryote_fdr_T)

result_isa_prokaryote_fdr_T[result_isa_prokaryote_fdr_T$s.NT==0 &
                         result_isa_prokaryote_fdr_T$s.T==1 ,] -> T

result_isa_prokaryote_fdr_T[result_isa_prokaryote_fdr_T$s.NT==1 &
                         result_isa_prokaryote_fdr_T$s.T==0 ,] -> NT
isa_df_T <- rbind(NT, T)
dim(isa_df_T)
isa_df_T

#Only SCN
result_isa_prokaryote_fdr_S <- isa_prokaryote_fdr_S$sign[which(isa_prokaryote_fdr_S$sign$p.value <= 0.05), ]
head(result_isa_prokaryote_fdr_S)
dim(result_isa_prokaryote_fdr_S)

result_isa_prokaryote_fdr_T[result_isa_prokaryote_fdr_S$s.Yes==0 &
                         result_isa_prokaryote_fdr_S$s.No==1 ,] -> No_SCN

result_isa_prokaryote_fdr_T[result_isa_prokaryote_fdr_S$s.Yes==1 &
                         result_isa_prokaryote_fdr_S$s.No==0 ,] -> Yes_SCN

isa_df_S <- rbind(Yes_SCN, No_SCN)
#isa_df <- rbind(NT, T, Yes_SCN, No_SCN, C, R, S, RCRC, SCSC, C_SCSC, S_SCSC, C_RCRC_SCSC, R_S_SCSC,R_RCRC_S_SCSC)
dim(isa_df_S)
isa_df_S

#Only ROTATION
result_isa_prokaryote_fdr_R <- isa_prokaryote_fdr_R$sign[which(isa_prokaryote_fdr_R$sign$p.value <= 0.05), ]
head(result_isa_prokaryote_fdr_R)
dim(result_isa_prokaryote_fdr_R)

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==1 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==0 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> R

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==1 &
                         result_isa_prokaryote_fdr_R$s.SCSC==0 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> C

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==1 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==0 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> S

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==1 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> SCSC

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==0 &
                         result_isa_prokaryote_fdr_R$s.RCRC==1 ,] -> RCRC

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==1 &
                         result_isa_prokaryote_fdr_R$s.SCSC==1 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> C_SCSC

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==1 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==1 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> S_SCSC

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==0 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==1 &
                         result_isa_prokaryote_fdr_R$s.SCSC==1 &
                         result_isa_prokaryote_fdr_R$s.RCRC==1 ,] -> C_RCRC_SCSC

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==1 &
                         result_isa_prokaryote_fdr_R$s.S==1 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==1 &
                         result_isa_prokaryote_fdr_R$s.RCRC==0 ,] -> R_S_SCSC

result_isa_prokaryote_fdr_R[result_isa_prokaryote_fdr_R$s.R==1 &
                         result_isa_prokaryote_fdr_R$s.S==0 &
                         result_isa_prokaryote_fdr_R$s.C==0 &
                         result_isa_prokaryote_fdr_R$s.SCSC==1 &
                         result_isa_prokaryote_fdr_R$s.RCRC==1 ,] -> R_RCRC_S_SCSC

isa_df_R <- rbind(R, C, S, SCSC, RCRC, C_SCSC, S_SCSC, C_RCRC_SCSC, R_S_SCSC,R_RCRC_S_SCSC)
#isa_df <- rbind(NT, T, Yes_SCN, No_SCN, C, R, S, RCRC, SCSC, C_SCSC, S_SCSC, C_RCRC_SCSC, R_S_SCSC,R_RCRC_S_SCSC)
dim(isa_df_R)
isa_df_R

# phyloseq objects of ISA OTUs --------------------------------------------------------------------
physeq_fungi -> physeq_fungi_isa_R
physeq_fungi_isa_R_ab = transform_sample_counts(physeq_fungi_isa_R, function(x) 100*x/sum(x)) # transform to relativa abundances 
otu_table(physeq_fungi_isa_R) <- otu_table(physeq_fungi_isa_R_ab)[rownames(isa_df_R), ]
physeq_fungi_isa_R
sample_data(physeq_fungi_isa_R)

physeq_prokaryote -> physeq_prokaryote_isa_R
physeq_prokaryote_isa_R_ab = transform_sample_counts(physeq_prokaryote_isa_R, function(x) 100*x/sum(x)) # transform to relativa abundances 
otu_table(physeq_prokaryote_isa_R) <-otu_table(physeq_prokaryote_isa_R_ab)[rownames(isa_df_R), ]
physeq_prokaryote_isa_R
sample_data(physeq_prokaryote_isa_R)

# improving taxonomy for indicator OTUs by BLAST
#write.dna(refseq(physeq_fungi_isa_R), format="fasta", colsep="", file="physeq_fungi_isa_R.fasta")
#write.csv(tax_table(physeq_fungi_isa_R), file = "tax_table_physeq_fungi_isa_R.csv")
#tax_physeq_fungi_isa_R <- tax_table(physeq_fungi_isa_R)
tax_table(physeq_fungi_isa_R) <- tax_table(as.matrix(tax_physeq_fungi_isa_R))
tax_fungi_isa_R <- tax_table(physeq_fungi_isa_R)

tax_table(physeq_prokaryote_isa_R) <- tax_table(as.matrix(tax_physeq_prokaryote_isa_R))
tax_prokaryote_isa_R <- tax_table(physeq_prokaryote_isa_R)
#write.csv(sample_data(physeq_fungi_isa_R) , file = "sample_data_physeq_fungi_isa_R.csv")
#sample_data_physeq_fungi_isa_R <- read.csv("sample_data_physeq_fungi_isa_R.csv", header=T, row.names =1)
#sample_data(physeq_fungi_isa_R) <- sample_data(sample_data_physeq_fungi_isa_R)
#sample_data(physeq_fungi_isa_R)$Isa <- factor(sample_data(physeq_fungi_isa_R)$Isa,
                                                 #levels=c("CapS1","CapS2","CapS3","CapS4","CapS5",
                                                          "CapS6","CapS7","CapS8","CapS9","CapS10",
                                                          "StemS1","StemS2","StemS3","StemS4","StemS5",
                                                          "StemS6","StemS7","StemS8","StemS9","StemS10", 
                                                          "SoilS1","SoilS2","SoilS3","SoilS4","SoilS5",
                                                          "SoilS6","SoilS7","SoilS8","SoilS9","SoilS10"))
#####
levels(sample_data(physeq_fungi_isa_R)$Rotation)
sample_data(physeq_fungi_isa_R)
otu_table(physeq_fungi_isa_R)

isa_obj_R <- as.data.frame(otu_table(physeq_fungi_isa_R))
dim(isa_obj_R)
isa_obj_R

isa_obj_T <- as.data.frame(otu_table(physeq_fungi_isa_T))
dim(isa_obj_T)
isa_obj_T

####
levels(sample_data(physeq_prokaryote_isa_R)$Rotation)
sample_data(physeq_prokaryote_isa_R)
otu_table(physeq_prokaryote_isa_R)

isa_obj_R <- as.data.frame(otu_table(physeq_prokaryote_isa_R))
dim(isa_obj_R)
isa_obj_R

isa_obj_T <- as.data.frame(otu_table(physeq_prokaryote_isa_T))
dim(isa_obj_T)
isa_obj_T

# Creating a data.frame for plotting the heatmap 
dentical(colnames(isa_obj_R), rownames(sample_data(physeq_fungi_isa_R)))
sample_data(physeq_fungi_isa_R)$Rotation
colnames(isa_obj_R) <- sample_data(physeq_fungi_isa_R)$Rotation
colnames(isa_df_R) <- c("C", "R", "RCRC","S", "SCSC", "Index", "Stat", "p.value")
identical(rownames(isa_obj_R), rownames(isa_df_R))

isa_obj_R <- cbind(isa_obj_R, isa_df_R)
isa_obj_R$readNo <- rowSums(otu_fungi_out[rownames(isa_df_R),])
isa_obj_R$relAb <- (isa_obj_R$readNo/sum(colSums(otu_fungi_out))) * 100
isa_obj_R$logAb <- log(isa_obj_R$readNo)
isa_obj_R$sqrtAb <- sqrt(isa_obj_R$readNo)
identical(rownames(tax_table(physeq_fungi_isa_R)), rownames(isa_obj_R))
isa_obj_R$taxa <- paste(rownames(isa_obj_R), as.data.frame(as.matrix(tax_table(physeq_fungi_isa_R)))$Rotation, sep = " ")
isa_obj_R$Rotation <- as.data.frame(as.matrix(tax_table(physeq_fungi_isa_R)))$Rotation
isa_obj_R

####
identical(colnames(isa_obj_R), rownames(sample_data(physeq_prokaryote_isa_R)))
sample_data(physeq_prokaryote_isa_R)$Rotation
colnames(isa_obj_R) <- sample_data(physeq_prokaryote_isa_R)$Rotation
colnames(isa_df_R) <- c("C", "R", "RCRC","S", "SCSC", "Index", "Stat", "p.value")
identical(rownames(isa_obj_R), rownames(isa_df_R))

isa_obj_R <- cbind(isa_obj_R, isa_df_R)
isa_obj_R$readNo <- rowSums(otu_prokaryote_out[rownames(isa_df_R),])
isa_obj_R$relAb <- (isa_obj_R$readNo/sum(colSums(otu_prokaryote_out))) * 100
isa_obj_R$logAb <- log(isa_obj_R$readNo)
isa_obj_R$sqrtAb <- sqrt(isa_obj_R$readNo)
identical(rownames(tax_table(physeq_prokaryote_isa_R)), rownames(isa_obj_R))
isa_obj_R$taxa <- paste(rownames(isa_obj_R), as.data.frame(as.matrix(tax_table(physeq_prokaryote_isa_R)))$Rotation, sep = " ")
isa_obj_R$Rotation <- as.data.frame(as.matrix(tax_table(physeq_prokaryote_isa_R)))$Rotation
isa_obj_R

# >>> HEATMAP -----------------------------------------------------------------------
library("ComplexHeatmap")
library("circlize")
library("seriation")

#BiocManager::install("ComplexHeatmap")
ex <- as.matrix(sqrt(isa_obj_R[,1:85]*10))
colnames(ex) <- factor(ex, levels = c("S", ""))



ht1 = Heatmap(ex, col = colorRamp2(c(0, 5), c("white","red")), 
              cluster_rows = TRUE, cluster_columns = FALSE, name = "Abundance",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = TRUE, row_labels = isa_obj_R$taxa)
ht1
#o = seriate(max(ht1) - ht1, method = "BEA_TSP")
#Heatmap(max(ht1) - ht1, name = "ht1", 
       # row_order = get_order(o, 1), column_order = get_order(o, 2))

ht2 = Heatmap(as.matrix(isa_obj_R[,31:33]), col = structure(c("red","white"), names = c("1","0")),
              cluster_rows = TRUE, cluster_columns = TRUE, name = "Group",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = TRUE, row_labels = tax_fungi_isa_R$Genus)

ha_bar = HeatmapAnnotation("Abundance (%)" = anno_barplot(isa_obj_R$sqrtAb/30, axis = TRUE, width = unit(2, "cm")), 
                           which = "row", show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 8),
                           annotation_name_rot = c(0)) # annotation_width = unit(20, "cm"), annotation_height = unit(10, "cm"),

# *** FIGURE 3 - HEATMAP -------------------------------------------------------------------------
ht1 + ha_bar + ht2 




# Creating a data.frame for plotting the heatmap 
dentical(colnames(isa_obj_T), rownames(sample_data(physeq_fungi_isa_T)))
sample_data(physeq_fungi_isa_T)$Rotation
colnames(isa_obj_T) <- sample_data(physeq_fungi_isa_T)$Rotation
colnames(isa_df_T) <- c("N", "NT", "Index", "Stat", "p.value")
identical(rownames(isa_obj_T), rownames(isa_df_T))

isa_obj_T <- cbind(isa_obj_T, isa_df_T)
isa_obj_T$readNo <- rowSums(otu_fungi_out[rownames(isa_df_T),])
isa_obj_T$relAb <- (isa_obj_T$readNo/sum(colSums(otu_fungi_out))) * 100
isa_obj_T$logAb <- log(isa_obj_T$readNo)
isa_obj_T$sqrtAb <- sqrt(isa_obj_T$readNo)
identical(rownames(tax_table(physeq_fungi_isa_T)), rownames(isa_obj_T))
isa_obj_T$taxa <- paste(rownames(isa_obj_T), as.data.frame(as.matrix(tax_table(physeq_fungi_isa_T)))$Rotation, sep = " ")
isa_obj_T$Rotation <- as.data.frame(as.matrix(tax_table(physeq_fungi_isa_T)))$Rotation
isa_obj_T

# >>> HEATMAP -----------------------------------------------------------------------
library("ComplexHeatmap")
library("circlize")
library("seriation")

#BiocManager::install("ComplexHeatmap")


ht1T= Heatmap(as.matrix(sqrt(isa_obj_T[,1:85]*10)), col = colorRamp2(c(0, 5), c("white","red")), 
              cluster_rows = TRUE, cluster_columns = FALSE, name = "Abundance",
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = TRUE, row_labels = isa_obj_T$taxa)
ht1T

# >>> COMMUNITY COMPOSITION **---------------------------------------------------------------------

# fungal composition -----------------------------------------------------------------------------
physeq_fungi_phylum_R <- tax_glom(physeq_fungi_isa_R, "Phylum")
otu_fungi_abund_phy_R <- taxa_sums(physeq_fungi_phylum_R)/sum(taxa_sums(physeq_fungi_phylum_R))*100
tax_abund_fungi_phy_R <- as(tax_table(physeq_fungi_phylum_R), "matrix")
tax_abund_fungi_phy_R <- as.data.frame(tax_abund_fungi_phy_R)
tax_abund_fungi_phy_R <- tax_abund_fungi_phy_R[c(2)]
tax_abund_fungi_phy_R$abundance <- as.vector(otu_fungi_abund_phy_R)
tax_abund_fungi_phy_R <- tax_abund_fungi_phy_R[order(tax_abund_fungi_phy_R$abundance, decreasing = TRUE),] 
tax_abund_fungi_phy_R
ggplot(tax_abund_fungi_phy_R, aes(x=reorder(Phylum, -abundance), y=abundance)) + 
  geom_bar(stat="identity") +
  labs(title="Soil", x="Taxa", y="Relative\nabundance %") +
  theme(axis.text.x = element_text(size=10, angle = -90, hjust = 0, vjust = 0.5))



physeq_fungi_order_R = tax_glom(physeq_fungi_isa_R, "Order")
otu_fungi_abund_ord_R = taxa_sums(physeq_fungi_order_R)/sum(taxa_sums(physeq_fungi_order_R))*100
tax_abund_fungi_ord_R <- as(tax_table(physeq_fungi_order_R), "matrix")
tax_abund_fungi_ord_R <- as.data.frame(tax_abund_fungi_ord_R)
tax_abund_fungi_ord_R <- tax_abund_fungi_ord_R[c(1:4)]
tax_abund_fungi_ord_R$abundance <- as.vector(otu_fungi_abund_ord_R)
tax_abund_fungi_ord_R <- tax_abund_fungi_ord_R[order(tax_abund_fungi_ord_R$abundance, decreasing = TRUE),] 
tax_abund_fungi_ord_R
ggplot(tax_abund_fungi_ord_R, aes(x=reorder(Order, -abundance), y=abundance)) + 
  geom_bar(stat="identity") +
  labs(title="Rotation", x="Taxa", y="Relative\nabundance %") +
  theme(axis.text.x = element_text(size=10, angle = -90, hjust = 0, vjust = 0.5))

physeq_fungi_genus_R = tax_glom(physeq_fungi_isa_R, "Genus")
otu_fungi_abund_gen_R = taxa_sums(physeq_fungi_genus_R)/sum(taxa_sums(physeq_fungi_genus_R))*100
tax_abund_fungi_gen_R <- as(tax_table(physeq_fungi_genus_R), "matrix")
tax_abund_fungi_gen_R <- as.data.frame(tax_abund_fungi_gen_R)
tax_abund_fungi_gen_R <- tax_abund_fungi_gen_R[c(1:6)]
tax_abund_fungi_gen_R$abundance <- as.vector(otu_fungi_abund_gen_R)
tax_abund_fungi_gen_R <- tax_abund_fungi_gen_R[order(tax_abund_fungi_gen_R$abundance, decreasing = TRUE),] 
tax_abund_fungi_gen_R
ggplot(tax_abund_fungi_gen_R, aes(x=reorder(Order, -abundance), y=abundance)) + 
  geom_bar(stat="identity") +
  labs(title="Rotation", x="Taxa", y="Relative\nabundance %") +
  theme(axis.text.x = element_text(size=10, angle = -90, hjust = 0, vjust = 0.5))

# extracting abundances for different stages
physeq_fungi_Young <- subset_samples(physeq_fungi, Stage%in%c("Young"))
otu_table(physeq_fungi_Young) <- otu_table(physeq_fungi_Young)[which(rowSums(otu_table(physeq_fungi_Young)) >= 1),] 
physeq_fungi_Young

physeq_fungi_Young_phylum = tax_glom(physeq_fungi_Young, "Phylum")
otu_fungi_abund_Young_phy = taxa_sums(physeq_fungi_Young_phylum)/sum(taxa_sums(physeq_fungi_Young_phylum))*100
tax_abund_fungi_Young_phy <- as(tax_table(physeq_fungi_Young_phylum), "matrix")
tax_abund_fungi_Young_phy <- as.data.frame(tax_abund_fungi_Young_phy)
tax_abund_fungi_Young_phy <- tax_abund_fungi_Young_phy[c(2)]
tax_abund_fungi_Young_phy$abundance <- as.vector(otu_fungi_abund_Young_phy)
tax_abund_fungi_Young_phy <- tax_abund_fungi_Young_phy[order(tax_abund_fungi_Young_phy$abundance, decreasing = TRUE),] 
tax_abund_fungi_Young_phy

physeq_fungi_Young_genus = tax_glom(physeq_fungi_Young, "Genus")
otu_fungi_abund_Young_gen = taxa_sums(physeq_fungi_Young_genus)/sum(taxa_sums(physeq_fungi_Young_genus))*100
tax_abund_fungi_Young_gen <- as(tax_table(physeq_fungi_Young_genus), "matrix")
tax_abund_fungi_Young_gen <- as.data.frame(tax_abund_fungi_Young_gen)
tax_abund_fungi_Young_gen <- tax_abund_fungi_Young_gen[c(1:6)]
tax_abund_fungi_Young_gen$abundance <- as.vector(otu_fungi_abund_Young_gen)
tax_abund_fungi_Young_gen <- tax_abund_fungi_Young_gen[order(tax_abund_fungi_Young_gen$abundance, decreasing = TRUE),] 
tax_abund_fungi_Young_gen[1:50, ]

physeq_fungi_Mature <- subset_samples(physeq_fungi, Stage%in%c("Mature"))
otu_table(physeq_fungi_Mature) <- otu_table(physeq_fungi_Mature)[which(rowSums(otu_table(physeq_fungi_Mature)) >= 1),] 
physeq_fungi_Mature

physeq_fungi_Mature_phylum = tax_glom(physeq_fungi_Mature, "Phylum")
otu_fungi_abund_Mature_phy = taxa_sums(physeq_fungi_Mature_phylum)/sum(taxa_sums(physeq_fungi_Mature_phylum))*100
tax_abund_fungi_Mature_phy <- as(tax_table(physeq_fungi_Mature_phylum), "matrix")
tax_abund_fungi_Mature_phy <- as.data.frame(tax_abund_fungi_Mature_phy)
tax_abund_fungi_Mature_phy <- tax_abund_fungi_Mature_phy[c(2)]
tax_abund_fungi_Mature_phy$abundance <- as.vector(otu_fungi_abund_Mature_phy)
tax_abund_fungi_Mature_phy <- tax_abund_fungi_Mature_phy[order(tax_abund_fungi_Mature_phy$abundance, decreasing = TRUE),] 
tax_abund_fungi_Mature_phy

physeq_fungi_Mature_genus = tax_glom(physeq_fungi_Mature, "Genus")
otu_fungi_abund_Mature_gen = taxa_sums(physeq_fungi_Mature_genus)/sum(taxa_sums(physeq_fungi_Mature_genus))*100
tax_abund_fungi_Mature_gen <- as(tax_table(physeq_fungi_Mature_genus), "matrix")
tax_abund_fungi_Mature_gen <- as.data.frame(tax_abund_fungi_Mature_gen)
tax_abund_fungi_Mature_gen <- tax_abund_fungi_Mature_gen[c(1:6)]
tax_abund_fungi_Mature_gen$abundance <- as.vector(otu_fungi_abund_Mature_gen)
tax_abund_fungi_Mature_gen <- tax_abund_fungi_Mature_gen[order(tax_abund_fungi_Mature_gen$abundance, decreasing = TRUE),] 
tax_abund_fungi_Mature_gen[1:50, ]

# prokaryotic composition ------------------------------------------------------------------------
#physeq_prokaryote_phylum = tax_glom(physeq_prokaryote, "Phylum")
#otu_prokaryote_abund_phy = taxa_sums(physeq_prokaryote_phylum)/sum(taxa_sums(physeq_prokaryote_phylum))*100
#tax_abund_prokaryote_phy <- as(tax_table(physeq_prokaryote_phylum), "matrix")
#tax_abund_prokaryote_phy <- as.data.frame(tax_abund_prokaryote_phy)
#tax_abund_prokaryote_phy <- tax_abund_prokaryote_phy[c(2)]
#tax_abund_prokaryote_phy$abundance <- as.vector(otu_prokaryote_abund_phy)
#tax_abund_prokaryote_phy <- tax_abund_prokaryote_phy[order(tax_abund_prokaryote_phy$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_phy

#physeq_prokaryote_class = tax_glom(physeq_prokaryote, "Class")
#otu_prokaryote_abund_class = taxa_sums(physeq_prokaryote_class)/sum(taxa_sums(physeq_prokaryote_class))*100
#tax_prokaryote_abund_class <- as(tax_table(physeq_prokaryote_class), "matrix")
#tax_prokaryote_abund_class <- as.data.frame(tax_prokaryote_abund_class)
#tax_prokaryote_abund_class <- tax_prokaryote_abund_class[c(1:3)]
#tax_prokaryote_abund_class$abundance <- as.vector(otu_prokaryote_abund_class)
#tax_prokaryote_abund_class <- tax_prokaryote_abund_class[order(tax_prokaryote_abund_class$abundance, decreasing = TRUE),] 
#tax_prokaryote_abund_class[1:50, ]

#physeq_prokaryote_order = tax_glom(physeq_prokaryote, "Order")
#otu_prokaryote_abund_ord = taxa_sums(physeq_prokaryote_order)/sum(taxa_sums(physeq_prokaryote_order))*100
#tax_prokaryote_abund_ord <- as(tax_table(physeq_prokaryote_order), "matrix")
#tax_prokaryote_abund_ord <- as.data.frame(tax_prokaryote_abund_ord)
#tax_prokaryote_abund_ord <- tax_prokaryote_abund_ord[c(1:4)]
#tax_prokaryote_abund_ord$abundance <- as.vector(otu_prokaryote_abund_ord)
#tax_prokaryote_abund_ord <- tax_prokaryote_abund_ord[order(tax_prokaryote_abund_ord$abundance, decreasing = TRUE),] 
#tax_prokaryote_abund_ord[1:50, ]

#physeq_prokaryote_genus = tax_glom(physeq_prokaryote, "Genus")
#otu_prokaryote_abund_gen = taxa_sums(physeq_prokaryote_genus)/sum(taxa_sums(physeq_prokaryote_genus))*100
#tax_abund_prokaryote_gen <- as(tax_table(physeq_prokaryote_genus), "matrix")
#tax_abund_prokaryote_gen <- as.data.frame(tax_abund_prokaryote_gen)
#tax_abund_prokaryote_gen <- tax_abund_prokaryote_gen[c(1:6)]
#tax_abund_prokaryote_gen$abundance <- as.vector(otu_prokaryote_abund_gen)
#tax_abund_prokaryote_gen <- tax_abund_prokaryote_gen[order(tax_abund_prokaryote_gen$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_gen[1:50, ]

# extracting abundances for different origins
#physeq_prokaryote_Cap <- subset_samples(physeq_prokaryote, Origin%in%c("Cap"))
#otu_table(physeq_prokaryote_Cap) <- otu_table(physeq_prokaryote_Cap)[which(rowSums(otu_table(physeq_prokaryote_Cap)) >= 1),] 
#physeq_prokaryote_Cap

#physeq_prokaryote_Cap_phylum = tax_glom(physeq_prokaryote_Cap, "Phylum")
#otu_prokaryote_abund_Cap_phy = taxa_sums(physeq_prokaryote_Cap_phylum)/sum(taxa_sums(physeq_prokaryote_Cap_phylum))*100
#tax_abund_prokaryote_Cap_phy <- as(tax_table(physeq_prokaryote_Cap_phylum), "matrix")
#tax_abund_prokaryote_Cap_phy <- as.data.frame(tax_abund_prokaryote_Cap_phy)
#tax_abund_prokaryote_Cap_phy <- tax_abund_prokaryote_Cap_phy[c(2)]
#tax_abund_prokaryote_Cap_phy$abundance <- as.vector(otu_prokaryote_abund_Cap_phy)
#tax_abund_prokaryote_Cap_phy <- tax_abund_prokaryote_Cap_phy[order(tax_abund_prokaryote_Cap_phy$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_Cap_phy

#physeq_prokaryote_Cap_genus = tax_glom(physeq_prokaryote_Cap, "Genus")
#otu_prokaryote_abund_Cap_gen = taxa_sums(physeq_prokaryote_Cap_genus)/sum(taxa_sums(physeq_prokaryote_Cap_genus))*100
#tax_abund_prokaryote_Cap_gen <- as(tax_table(physeq_prokaryote_Cap_genus), "matrix")
#tax_abund_prokaryote_Cap_gen <- as.data.frame(tax_abund_prokaryote_Cap_gen)
#tax_abund_prokaryote_Cap_gen <- tax_abund_prokaryote_Cap_gen[c(1:6)]
#tax_abund_prokaryote_Cap_gen$abundance <- as.vector(otu_prokaryote_abund_Cap_gen)
#tax_abund_prokaryote_Cap_gen <- tax_abund_prokaryote_Cap_gen[order(tax_abund_prokaryote_Cap_gen$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_Cap_gen[1:50, ]

#physeq_prokaryote_Stem <- subset_samples(physeq_prokaryote, Origin%in%c("Stem"))
#otu_table(physeq_prokaryote_Stem) <- otu_table(physeq_prokaryote_Stem)[which(rowSums(otu_table(physeq_prokaryote_Stem)) >= 1),] 
#physeq_prokaryote_Stem

#physeq_prokaryote_Stem_phylum = tax_glom(physeq_prokaryote_Stem, "Phylum")
#otu_prokaryote_abund_Stem_phy = taxa_sums(physeq_prokaryote_Stem_phylum)/sum(taxa_sums(physeq_prokaryote_Stem_phylum))*100
#tax_abund_prokaryote_Stem_phy <- as(tax_table(physeq_prokaryote_Stem_phylum), "matrix")
#tax_abund_prokaryote_Stem_phy <- as.data.frame(tax_abund_prokaryote_Stem_phy)
#tax_abund_prokaryote_Stem_phy <- tax_abund_prokaryote_Stem_phy[c(2)]
#tax_abund_prokaryote_Stem_phy$abundance <- as.vector(otu_prokaryote_abund_Stem_phy)
#tax_abund_prokaryote_Stem_phy <- tax_abund_prokaryote_Stem_phy[order(tax_abund_prokaryote_Stem_phy$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_Stem_phy

#physeq_prokaryote_Stem_genus = tax_glom(physeq_prokaryote_Stem, "Genus")
#otu_prokaryote_abund_Stem_gen = taxa_sums(physeq_prokaryote_Stem_genus)/sum(taxa_sums(physeq_prokaryote_Stem_genus))*100
#tax_abund_prokaryote_Stem_gen <- as(tax_table(physeq_prokaryote_Stem_genus), "matrix")
#tax_abund_prokaryote_Stem_gen <- as.data.frame(tax_abund_prokaryote_Stem_gen)
#tax_abund_prokaryote_Stem_gen <- tax_abund_prokaryote_Stem_gen[c(1:6)]
#tax_abund_prokaryote_Stem_gen$abundance <- as.vector(otu_prokaryote_abund_Stem_gen)
#tax_abund_prokaryote_Stem_gen <- tax_abund_prokaryote_Stem_gen[order(tax_abund_prokaryote_Stem_gen$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_Stem_gen[1:50, ]

#physeq_prokaryote_Soil <- subset_samples(physeq_prokaryote, Origin%in%c("Soil"))
#otu_table(physeq_prokaryote_Soil) <- otu_table(physeq_prokaryote_Soil)[which(rowSums(otu_table(physeq_prokaryote_Soil)) >= 1),] 
#physeq_prokaryote_Soil

#physeq_prokaryote_Soil_phylum = tax_glom(physeq_prokaryote_Soil, "Phylum")
#otu_prokaryote_abund_Soil_phy = taxa_sums(physeq_prokaryote_Soil_phylum)/sum(taxa_sums(physeq_prokaryote_Soil_phylum))*100
#tax_abund_prokaryote_Soil_phy <- as(tax_table(physeq_prokaryote_Soil_phylum), "matrix")
#tax_abund_prokaryote_Soil_phy <- as.data.frame(tax_abund_prokaryote_Soil_phy)
#tax_abund_prokaryote_Soil_phy <- tax_abund_prokaryote_Soil_phy[c(2)]
#tax_abund_prokaryote_Soil_phy$abundance <- as.vector(otu_prokaryote_abund_Soil_phy)
#tax_abund_prokaryote_Soil_phy <- tax_abund_prokaryote_Soil_phy[order(tax_abund_prokaryote_Soil_phy$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_Soil_phy

#physeq_prokaryote_Soil_genus = tax_glom(physeq_prokaryote_Soil, "Genus")
#otu_prokaryote_abund_Soil_gen = taxa_sums(physeq_prokaryote_Soil_genus)/sum(taxa_sums(physeq_prokaryote_Soil_genus))*100
#tax_abund_prokaryote_Soil_gen <- as(tax_table(physeq_prokaryote_Soil_genus), "matrix")
#tax_abund_prokaryote_Soil_gen <- as.data.frame(tax_abund_prokaryote_Soil_gen)
#tax_abund_prokaryote_Soil_gen <- tax_abund_prokaryote_Soil_gen[c(1:6)]
#tax_abund_prokaryote_Soil_gen$abundance <- as.vector(otu_prokaryote_abund_Soil_gen)
#tax_abund_prokaryote_Soil_gen <- tax_abund_prokaryote_Soil_gen[order(tax_abund_prokaryote_Soil_gen$abundance, decreasing = TRUE),] 
#tax_abund_prokaryote_Soil_gen[1:50, ]

# looking for specific genera
tax_abund_prokaryote_Soil_gen[tax_abund_prokaryote_Soil_gen$Genus=="Pedobacter",]

# plottign bars ----------------------------------------------------------------------------------
library("scales")
library("grid")
library("reshape2")

# melt to long format (for ggploting) ------------------------------------------------------------
# and  prune out phyla below 2% in each sample
otu_fungi_phylum <- physeq_fungi %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

otu_fungi_phylum
head(otu_fungi_phylum)
unique(otu_fungi_phylum$Family)

# Set colors for plotting
phylum_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5",
                   "#5E738F","#D1A33D", "#8A7C64", "#599861")

# *** FIGURE 1A barplot FUNGI 
barplot_fungi = ggplot(otu_fungi_phylum, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") +
  labs(title="", x="Samples", y = "Relative Abundance") +
  scale_fill_discrete() +
  facet_grid(~Rotation, scales = "free_x", space="free_x") +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 10)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")

barplot_fungi


# composition Prokaryotic communities ------------------------------------------------------------
#sample_data(physeq_prokaryote)$Origin <- factor(sample_data(physeq_prokaryote)$Origin,
   #                                             levels=c("Cap","Stem","Soil"))

#otu_prokaryote_phylum <- physeq_prokaryote %>%
#  tax_glom(taxrank = "Phylum") %>%                     
#  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
#  psmelt() %>%                                         
#  filter(Abundance > 0.01) %>%                         
#  arrange(Phylum)                                     

#otu_prokaryote_phylum
#head(otu_prokaryote_phylum)
#unique(otu_prokaryote_phylum$Phylum)

# *** FIGURE 1B barplot prokaryote 
#barplot_prokaryote = ggplot(otu_prokaryote_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
#  geom_bar(stat = "identity") +
 # labs(title="", x="Samples", y = "Relative Abundance") +
 # scale_fill_manual(values = phylum_colors) +
 # facet_grid(~Origin, scales = "free_x", space="free_x") +
#  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
#  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
#  theme(strip.text.x = element_text(size = 8, face = "bold")) +
 # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
 # theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
 # theme(legend.position="right")

#barplot_prokaryote

# *** FIGURE 1 - BARPLOTS ------------------------------------------------------------------------

ggarrange(barplot_fungi,
          barplot_prokaryote,
          labels = c("A", "B"),
          widths = c(1.02,2),
          align = "none", 
          ncol = 2, nrow = 1,
          common.legend = FALSE,
          legend = c("right"))

# >>> VENN DIAGRAM -------------------------------------------------------------------------------
library("limma")
source("my_venn_diag.R")

physeq_fungi_R = merge_samples(physeq_fungi, "Rotation")
otu_fungi_R <- as.data.frame(t(otu_table(physeq_fungi_R)))
venn_counts_otu_fungi_R <- vennCounts(otu_fungi_R, include="both")
venn_counts_otu_fungi_R

physeq_prokaryote_R = merge_samples(physeq_prokaryote, "Rotation")
otu_prokaryote_R <- as.data.frame(t(otu_table(physeq_prokaryote_R)))
venn_counts_otu_prokaryote_R <- vennCounts(otu_prokaryote_R, include="both")
venn_counts_otu_prokaryote_R
#physeq_prokaryote_St = merge_samples(physeq_prokaryote, "Stage")
#otu_prokaryote_St <- as.data.frame(t(otu_table(physeq_prokaryote_St)))
#venn_counts_otu_prokaryote_St <- vennCounts(otu_prokaryote_St, include="both")
#venn_counts_otu_prokaryote_St

#physeq_prokaryote_Ts = merge_samples(physeq_prokaryote, "Origin")
#otu_prokaryote_Ts <- as.data.frame(t(otu_table(physeq_prokaryote_Ts)))
#venn_counts_otu_prokaryote_Ts <- vennCounts(otu_prokaryote_Ts, include="both")
#venn_counts_otu_prokaryote_Ts

# *** FIGURE 4 - venn diagrams -------------------------------------------------------------------

#layout(matrix(1:3, ncol=3))
#venn_Gian(venn_counts_otu_prokaryote_Ts,
          #cex=c(1),
         # circle.col =c("#bdbdbd", "#525252", "#000000"),
         # mar = c(1,1,1,1),
         # lwd = 2, main="A")

#venn_Gian(venn_counts_otu_prokaryote_St,
          #cex=c(1),
          #circle.col =c("red", "grey"),
          #mar = c(1,1,1,1),
          #lwd = 2, main="B")

venn_Gian(venn_counts_otu_fungi_R,
          cex=c(1),
          circle.col =c("red", "grey","blue","green","purple"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")

venn_Gian(venn_counts_otu_prokaryote_R,
          cex=c(1),
          circle.col =c("red", "grey","blue","green","purple"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")
#dev.off()

#Tillage
physeq_fungi_T = merge_samples(physeq_fungi, "Tillage")
otu_fungi_T <- as.data.frame(t(otu_table(physeq_fungi_T)))
venn_counts_otu_fungi_T <- vennCounts(otu_fungi_T, include="both")
venn_counts_otu_fungi_T

physeq_prokaryote_T = merge_samples(physeq_prokaryote, "Tillage")
otu_prokaryote_T <- as.data.frame(t(otu_table(physeq_prokaryote_T)))
venn_counts_otu_prokaryote_T <- vennCounts(otu_prokaryote_T, include="both")
venn_counts_otu_prokaryote_T

#physeq_prokaryote_St = merge_samples(physeq_prokaryote, "Stage")
#otu_prokaryote_St <- as.data.frame(t(otu_table(physeq_prokaryote_St)))
#venn_counts_otu_prokaryote_St <- vennCounts(otu_prokaryote_St, include="both")
#venn_counts_otu_prokaryote_St

#physeq_prokaryote_Ts = merge_samples(physeq_prokaryote, "Origin")
#otu_prokaryote_Ts <- as.data.frame(t(otu_table(physeq_prokaryote_Ts)))
#venn_counts_otu_prokaryote_Ts <- vennCounts(otu_prokaryote_Ts, include="both")
#venn_counts_otu_prokaryote_Ts

# *** FIGURE 4 - venn diagrams -------------------------------------------------------------------

#layout(matrix(1:3, ncol=3))
#venn_Gian(venn_counts_otu_prokaryote_Ts,
         # cex=c(1),
       #   circle.col =c("#bdbdbd", "#525252", "#000000"),
         # mar = c(1,1,1,1),
         # lwd = 2, main="A")

#venn_Gian(venn_counts_otu_prokaryote_St,
#cex=c(1),
#circle.col =c("red", "grey"),
#mar = c(1,1,1,1),
#lwd = 2, main="B")

venn_Gian(venn_counts_otu_fungi_T,
          cex=c(1),
          circle.col =c("red", "grey","blue","green","purple"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")

venn_Gian(venn_counts_otu_prokaryote_T,
          cex=c(1),
          circle.col =c("red", "grey","blue","green","purple"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")
#dev.off()

#SCN
physeq_fungi_S = merge_samples(physeq_fungi, "SCN")
otu_fungi_S <- as.data.frame(t(otu_table(physeq_fungi_S)))
venn_counts_otu_fungi_S <- vennCounts(otu_fungi_S, include="both")
venn_counts_otu_fungi_S

physeq_prokaryote_S = merge_samples(physeq_prokaryote, "SCN")
otu_prokaryote_S <- as.data.frame(t(otu_table(physeq_prokaryote_S)))
venn_counts_otu_prokaryote_S <- vennCounts(otu_prokaryote_S, include="both")
venn_counts_otu_prokaryote_S

#physeq_prokaryote_St = merge_samples(physeq_prokaryote, "Stage")
#otu_prokaryote_St <- as.data.frame(t(otu_table(physeq_prokaryote_St)))
#venn_counts_otu_prokaryote_St <- vennCounts(otu_prokaryote_St, include="both")
#venn_counts_otu_prokaryote_St

#physeq_prokaryote_Ts = merge_samples(physeq_prokaryote, "Origin")
#otu_prokaryote_Ts <- as.data.frame(t(otu_table(physeq_prokaryote_Ts)))
#venn_counts_otu_prokaryote_Ts <- vennCounts(otu_prokaryote_Ts, include="both")
#venn_counts_otu_prokaryote_Ts

# *** FIGURE 4 - venn diagrams -------------------------------------------------------------------

#layout(matrix(1:3, ncol=3))
#venn_Gian(venn_counts_otu_prokaryote_Ts,
# cex=c(1),
#   circle.col =c("#bdbdbd", "#525252", "#000000"),
# mar = c(1,1,1,1),
# lwd = 2, main="A")

#venn_Gian(venn_counts_otu_prokaryote_St,
#cex=c(1),
#circle.col =c("red", "grey"),
#mar = c(1,1,1,1),
#lwd = 2, main="B")

venn_Gian(venn_counts_otu_fungi_S,
          cex=c(1),
          circle.col =c("red", "grey","blue","green","purple"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")

venn_Gian(venn_counts_otu_prokaryote_S,
          cex=c(1),
          circle.col =c("red", "grey","blue","green","purple"),
          mar = c(1,1,1,1),
          lwd = 2, main="C")
#dev.off()

# >>> MICROBIAL NETWORK - prokaryotic communitie -----------------------------------------------------

# Note: This is just an example of the code used to plot the 
# networks included in the manuscript. For this reason, the network
# figure (the plot) will not be the exact same of the one included 
# in the manuscript.

# extarcting taxa present in >= 17 samples. This value is choosen 
# according to: 1) group of samples present in the variable with the most 
# sample per treatment, 2) by visually estimating the best level of 
# sparsity that will make clear to distinguish  important network
# properties, and 3) by an optimal network stability.

physeq_fungi -> physeq_fungi_filt
cntNonZero <- apply(as.data.frame(otu_table(physeq_fungi)), 1, function(x) sum(x > 20))
otu_table(physeq_fungi_filt) <- otu_table(physeq_fungi_filt)[which(cntNonZero >=20),]
rowSums(otu_table(physeq_fungi_filt))
physeq_fungi_filt

# Network analysis -------------------------------------------------------------------------------
#library(devtools)
#install_github("zdk123/SpiecEasi")
library("SpiecEasi"); packageVersion("SpiecEasi")
library("igraph"); packageVersion("igraph")
library("huge"); packageVersion("huge")
library("qgraph"); packageVersion("qgraph")
library("MCL")

set.seed(2020)

spiec_prok <- spiec.easi(biom_16S_uparse_filt, 
                         method="mb",
                         lambda.min.ratio=1e-2, 
                         nlambda=50, 
                         sel.criterion ="stars", 
                         pulsar.select = TRUE,
                         pulsar.params=list(rep.num=90))

spiec_prok
getStability(spiec_prok)

??spiec.easi

# creating the network object 
plot_spiec_prok <- adj2igraph(getRefit(spiec_prok),
                              vertex.attr=list(name=taxa_names(biom_16S_uparse_filt)))

plot_spiec_prok


# PLOT function ----------------------------------------------------------------------------------
# modified from:
# Beiko, R. G., Hsiao, W., and Parkinson, J.  
# Microbiome Analysis: Methods and Protocols.
# Humana Press, 2018.

plot.net.cls <- function(net, scores, cls, art_point, physeq_obj) {
  # Get size of clusters to find isolated nodes.
  cls_sizes <- sapply(groups(cls), length)
  # Randomly choosing node colors. Users can provide their own vector of colors.
  colors <- sample(colours(), length(cls))
  # Nodes in clusters will be color coded. Isolated nodes will be white.
  V(net)$color <- sapply(membership(cls),
                         function(x) {ifelse(cls_sizes[x]>1,
                                             colors[x], "white")})
  # Convert node label from names to numerical IDs.
  node.names <- V(net)$name
  #col_ids <- seq(1, length(node.names))
  #V(net)$name <- col_ids
  # name nodes with taxonomic names 
  #taxa_physeq_obj <- as.data.frame(as.matrix(tax_table(physeq_obj)))
  #V(net)$name <- taxa_physeq_obj$Phylum
  # To draw a halo around articulation points.
  art_point <- lapply(names(art_point), function(x) x)
  marks <- lapply(1:length(art_point), function(x) which(node.names == art_point[[x]]))
  # vertex size 
  sum_list <- taxa_sums(physeq_obj)
  vsize = log(sum_list)/2
  # set size of vertex proportional to clr-mean
  #vsize <- rowMeans(clr(t(otu_table(physeq_obj)), 1)) + 7
  # Customized layout to avoid nodes overlapping.
  #e <- get.edgelist(net)
  #class(e) <- "numeric"
  #l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(net),
  #                                         area=8*(vcount(net)^2),
  #                                        repulse.rad=(vcount(net)^3.1))
  # Main plot function.
  plot.igraph(net, 
              vertex.size = vsize, 
              vertex.label.cex=0.4,
              vertex.label.color = "black",
              mark.border="black",
              mark.groups = marks,
              mark.col = "white",
              mark.expand = 10,
              mark.shape = 1,
              layout = NULL)
}

# *** FIGURE 5 - the network ---------------------------------------------------------------------
plot.net.cls(plot_spiec_prok, 
             hub_score(plot_spiec_prok)$vector, 
             walktrap.community(plot_spiec_prok), 
             articulation.points(plot_spiec_prok), 
             physeq_fungi_filt)

# >>> NETOWRK TOPOLOGY INDEXES ------------------------------------------------------------------
cent_res <- igraph::centrality(plot_spiec_prok, all.shortest.paths = TRUE)
cent_res$OutDegree
cent_res$InDegree
cent_res$Closeness
cent_res$Betweenness

# Degree 
igraph::degree(plot_spiec_prok, mode="all")

# Degree distribution 
deg.dist1 <- degree_distribution(plot_spiec_prok, mode = "all")
deg.dist1
plot(deg.dist1, type='b', ylab="Frequency", xlab="Degree", main="Degree Distribution")

# Transitivity or clustering coefficient
clustering_coeff <- transitivity(plot_spiec_prok, type = "global")
clustering_coeff


# _______________ IMPORTANT NOTE ___________________-----------------

# To relabel the Origin factor with the correct scientific terminology 
# we used the code below and then we re-sun the R scripts 
sample_data(physeq_prokaryote_mSeq)$Origin <- c(rep("Pileus", 10),rep("Soil",10), rep("Stipe",10)) 
sample_data(physeq_prokaryote_mSeq)$Origin <- factor(sample_data(physeq_prokaryote_mSeq)$Origin,levels=c("Pileus","Stipe","Soil"))
