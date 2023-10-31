# Charli Davies

# R code and analysis for: Immunogenetic variation shapes the gut microbiome in a natural vertebrate population

# last modified 09/11/2021

# using R version 3.6.1

# Need to first put all files from dryad into the same folder and set this as you're working directory

# clear R of all objects
rm(list=ls())

#load libraries####
library(phyloseq)
library(ggplot2)
library(ape)
library(lme4)
library(tidyverse)
library(btools)
library(FSA)
library(gridExtra)
library(grid)

# set working directory####

setwd("~/Decontam")

# import data into R and prepare for phyloseq####

# import OTU table
otu<-read.table("otu_table-inccontam.txt", header=T)
head(otu)
View(otu)

#import taxonomy

taxon<-read.csv("taxonomy.tsv", sep ='\t', header = T)
head(taxon)
View(taxon)

# merge the two files - this will remove any columns present in taxonomy but not otu table - ie chloroplasts etc
merged<-merge(otu, taxon, by= "OTUID")
head(merged)
View(merged)

# output the merged file - then need to seperate into two files again - one for taxonomy and one for otu table - for taxonomy file need to seperate text > columns to seperate out kingdom, phylum etc and save both files as .csv

write.table(merged, "merged_otu_table_contam.txt", sep="\t", col.names=TRUE, row.names=FALSE)


# once files changed can start with phyloseq



##################################### Import into phyloseq ########################################

# Import data ####
# Import otu table

otu<-read.csv("OTU_TABLE_CONTAM.csv", sep=",", row.names=1) 
otu_table<-as.matrix(otu)
# View(otu)

# IMPORT TAXONOMY

taxonomy<-read.csv("Taxonomy_silva_CONTAM.csv", sep=",", row.names=1)
taxonomy<-as.matrix(taxonomy)
# View(taxonomy)


# Import metadata
metadata<-read.csv("Sample_metadata_Daviesetal.csv", sep=",", row.names = 1)
View(metadata)

#Import rooted tree
phy_tree<-read_tree("tree.nwk")
is.rooted(phy_tree(phy_tree))

# need to import csv documents as phyloseq documents (tree already is a phyloseq document)
OTU<-otu_table(otu, taxa_are_rows = TRUE) 
TAX<-tax_table(taxonomy)
META<-sample_data(metadata)

#check names are consistent accross objects
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

# make sure files have same sample names
sample_names(OTU)
sample_names(META)

# merge into one phyloseq object
physeq<- phyloseq(OTU,TAX,META,phy_tree)
summary(physeq)
physeq


# check for otus that are not present in any sample - none 
any(sample_sums(physeq) == 0)

# plots out total reads per sample and distribution
readsumsdf = data.frame(nreads = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq),type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(physeq),TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# Subset dataset + reorder####

#remove sample with no reads and ASVs WITH less than 50 abundance
physeq <-prune_taxa(taxa_sums(physeq) > 50, physeq) 
physeq <-prune_samples(sample_sums(physeq)>=0, physeq)
physeq <- subset_samples(physeq, CleanedFeatureCount != "0")

#First only include negative or faecal samples
Contam <- subset_samples(physeq, SampleType4way%in%c("F","Extraction control"))
Control <- subset_samples(physeq, SampleType4way%in%c("Extraction control"))

# Merge genus
genus_Contam <- tax_glom(Contam, taxrank="Genus")
#extract taxonomy and OTU counts from merged genus data

#install 
library(decontam); packageVersion("decontam")

# make plot of reads
df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType4way)) + geom_point()

# identify contaminants by prevelance
sample_data(genus_Contam)$is.neg <- sample_data(genus_Contam)$SampleType4way == "Extraction control"
contamdf.prev <- isContaminant(genus_Contam, method="prevalence", neg="is.neg", batch = "RunID")
table(contamdf.prev$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(genus_Contam, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleType4way == "Extraction control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleType4way == "F", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Can also identify by frequency - here will use cleaned feature count as measure of frequency

contamdf.freq <- isContaminant(genus_Contam, method="frequency", conc="CleanedFeatureCount", batch = "RunID",neg="is.neg")
head(contamdf.freq)
table(contamdf.freq$contaminant)
View(contamdf.freq$contaminant$True)

plot_frequency(genus_Contam, taxa_names(genus_Contam)[c(1,2,3)], conc="CleanedFeatureCount") + 
  xlab("Number of reads")

# using both methods together

Both <- isContaminant(genus_Contam, conc="CleanedFeatureCount", neg="is.neg", detailed=TRUE, normalize=TRUE, method='both', batch = "RunID")
table(Both$contaminant)
view(Both$contaminant)

### remove contaminants and find out what they are
contaminants<- row.names(subset(Both, contaminant == TRUE))
only_contaminants<-prune_taxa(contaminants, Contam)
taxa.names<-taxa_names(only_contaminants)
View(tax_table(only_contaminants)@.Data)




ntaxa(Control)
sum(sample_sums(Control))

# Subset out just the negative controls. 


# Clean out taxa/SV columns that are no longer present.
ps_ex_controls <- prune_taxa(taxa_sums(Control) > 0, Control)

# How abundant are the SVs in the negative controls?	
plot(sort(taxa_sums(ps_ex_controls), TRUE), type="h")

# What's in them?	
plot_bar(ps_ex_controls, "UniqueNumber", fill = "Genus")

# Find taxa that are shared among a certain number of negative controls.  This might indicate if there is a common contaminant that might be coming from a common reagent, often water.
controls_shared_taxa = filter_taxa(ps_ex_controls, function(x) sum(x >= 1) == (6), TRUE)
plot_bar(controls_shared_taxa, "UniqueNumber", fill = "Genus")
