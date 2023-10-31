# Charli Davies

# R code and analysis for: Immunogenetic variation shapes the gut microbiome in a natural vertebrate population

# last modified 09/11/2021

# using R version 3.6.1

# Need to first put all files from dryad into the same folder and set this as you're working directory
# Note some files have been exported from R, sorted in excel for aesthetic purposes and then re-exported in R - see folder for changed files in github

# clear R of all objects
rm(list=ls())

#load libraries - also load libraries throughout script####
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

setwd("~")

# import data into R and prepare for phyloseq####

# import OTU table
otu<-read.table("otu_table.txt", header=T)
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

# output the merged file 

write.table(merged, "merged_otu_table.txt", sep="\t", col.names=TRUE, row.names=FALSE)

# then need to seperate into two files again - one for taxonomy and one for otu table - for taxonomy file need to seperate text > columns to seperate out kingdom, phylum etc and remove D_ and save files as OTU_TABLE.csv and Taxonomy_silva.csv

# once files changed can start with phyloseq


##################################### Import data into phyloseq ########################################

# Import data ####
# Import otu table

otu<-read.csv("OTU_TABLE.csv", sep=",", row.names=1) 
otu_table<-as.matrix(otu)
# View(otu)

# IMPORT TAXONOMY

taxonomy<-read.csv("Taxonomy_silva.csv", sep=",", row.names=1)
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

#Re-order age category variables
correct.order <- c("CH", "FL", "OFL","SA", "A", "NA")
sample_data(physeq)$AgeClass5 <- factor(sample_data(physeq)$AgeClass5,levels = correct.order)

# check for otus that are not present in any sample - none 
any(sample_sums(physeq) == 0)

# plots out total reads per sample and distribution
readsumsdf = data.frame(nreads = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq),type = "OTUs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(physeq),TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# Subset dataset + reorder####

#remove sample with no reads
physeq <-prune_samples(sample_sums(physeq)>=10, physeq)

#only faecal samples
Faecal <- subset_samples(physeq, SampleType%in%c("F"))

# only repeated samples - for measures of similarity
Repeat <- subset_samples(Faecal, SampleRepeat%in%c("Yes"))

#ONLY COUSIN faecal samples (also remove samples sequenced twice and samples from the same catch)
Island <- subset_samples(Faecal, Island%in%c("CN"))
Cousin <- subset_samples(Island, Keep_Catch_Bag_Tray%in%c("y"))

# only samples with read depth greater than 10000
Cousin_final <-prune_samples(sample_sums(Cousin)>=10000, Cousin)

#remove samples with no MHC data
Cousin_final_mhc <-subset_samples(Cousin_final, Both_MHC_present%in%c("Y"))

#remove chick samples
Cousin_final_mhc_nochick <-subset_samples(Cousin_final_mhc, Chick%in%c("NotChick"))

#remove duplicate birds samples
Cousin_final_mhc_nochick_nodup <-subset_samples(Cousin_final_mhc_nochick, Keep_DUPLICATE_BIRD%in%c("y"))

# only ASVs with read depth greater than 50
Cousin_final_cleaned <-prune_taxa(taxa_sums(Cousin_final_mhc_nochick_nodup)>50,Cousin_final_mhc_nochick_nodup)

# to rarefy samples to 10000 reads
set.seed(28367)
Cousin_final_cleaned_rare = rarefy_even_depth(Cousin_final_cleaned, sample.size = 10000)

# To calculate SEM
standard_error <- function(x) sd(x) / sqrt(length(x))

#How Many Taxa/sequences did We Lose in rarefaction?    
ntaxa(Cousin_final); ntaxa(Cousin_final_cleaned); ntaxa(Cousin_final_cleaned_rare) #34869/9628/ 9614
sum(sample_sums(Cousin_final)); sum(sample_sums(Cousin_final_cleaned)); sum(sample_sums(Cousin_final_cleaned_rare)) # 6562592/10680281/ 1950000
mean(sample_sums(Cousin_final)); mean(sample_sums(Cousin_final_cleaned)); mean(sample_sums(Cousin_final_cleaned_rare)) # 58941.6/54770.7/ 10000
standard_error(sample_sums(Cousin_final)); standard_error(sample_sums(Cousin_final_cleaned)); standard_error(sample_sums(Cousin_final_cleaned_rare)) #3997.7/4170.4/0

# Merge by taxa####
phylum_merge <- tax_glom(Cousin_final, taxrank="Phylum")
class_merge <- tax_glom(Cousin_final, taxrank="Class") 
order_merge <- tax_glom(Cousin_final, taxrank="order") 
family_merge <- tax_glom(Cousin_final, taxrank="Family") 
genus_merge <- tax_glom(Cousin_final, taxrank="Genus")
species_merge <- tax_glom(Cousin_final, taxrank="Species")

#to see how many of each taxa
ntaxa(phylum_merge)
ntaxa(class_merge)
ntaxa(order_merge)
ntaxa(family_merge)
ntaxa(genus_merge)
ntaxa(species_merge)
ntaxa(Cousin_final)


##################################### Fig S1 ####

# Sample completeness and alpha rarefaction plotting (species accumulation) using the package iNEXT 

head(otu_table(Cousin_final)) #check you can access OTU abundance table in phyloseq object
str(otu_table(Cousin_final))

abundances<-as(otu_table(Cousin_final), "matrix") #make the otu abundance table into a matrix (iNEXT only accepts a particular format as input)
abundances[1:2,1:2] #check the matrix

abundances2<-as.data.frame(abundances) #change to a dataframe
str(abundances2)

library(dplyr)
df2 <- mutate_all(abundances2, function(x) as.numeric(x)) #iNEXT only takes numeric values, so change all values in the dataframe to numeric values instead of integers.
str(df2)

#install.packages("iNEXT")
library(iNEXT)
library(ggplot2)

inext<-iNEXT(df2, q=0, datatype="abundance", endpoint=20000) #runs the inext function on your dataframe. 

#plot a rarefaction curve using species richness:
rarefaction1<- ggiNEXT(inext, type=1, se=TRUE, facet.var="none", grey=TRUE)+theme(legend.position = "none")+ xlab("Sequencing depth") + ylab("Observed ASVs")+theme(axis.title.x=element_text(size=10, color="black"),axis.text.x = element_text(size=10,color="black" ), axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))+
  theme(axis.title.y=element_text(color="black", size=10 , margin=margin(0,5,0,0)),axis.text.y=element_text(size=10,color="black" ), axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))
rarefaction<-rarefaction1+ geom_vline(xintercept=10000, alpha=0.5, linetype=2)

#plot a sample completeness curve:
completeness1<-ggiNEXT(inext, type=2, se=TRUE, facet.var="none", grey=TRUE)+scale_shape_manual(values=rep(20,32))+ theme(legend.position = "none") +xlab("Read count") +ylab("Sample completeness")+theme(axis.title.x=element_text(size=10, color="black"),axis.text.x = element_text(size=10,color="black" ), axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))+
  theme(axis.title.y=element_text(color="black", size=10 , margin=margin(0,5,0,0)),axis.text.y=element_text(size=10,color="black" ), axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))
completeness <-completeness1+ geom_vline(xintercept=10000, alpha=0.5, linetype=2)

# add in labels
completenessplot <- arrangeGrob(completeness, top = textGrob("a",  x = unit(0, "npc"), y   = unit(0.4, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=10)))
rarefactionplot <- arrangeGrob(rarefaction, top = textGrob("b",  x = unit(0, "npc"), y   = unit(0.4, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=10)))

#export graph
combo_plot<-grid.arrange(completenessplot,rarefactionplot, ncol = 2)
ggsave("FigS1.jpeg", combo_plot, width = 16, height =6 , units="cm", dpi = 300)


##################################### Figure 1 ########

#A:phylum graph####
#this gets down into merged columns to make tidier graph have additionally re-ordered in excel

phylum_merge <- tax_glom(Cousin_final, taxrank="Phylum") 
View(phylum_merge@tax_table@.Data)

# Find top 10 phyla by abundance
top10phyla <- sort(taxa_sums(phylum_merge),decreasing=TRUE)[1:10] %>% names()
rarephyla <- setdiff(taxa_names(phylum_merge),top10phyla)
rarephyla
top10phyla

# Merge rare phyla into "Low abundance" category
merged.top.phylum<-merge_taxa(phylum_merge,rarephyla,archetype="2480e81716ba1920a54162ee3edbcc99") %>%  prune_taxa(taxa_sums(.)>0,.)
ntaxa(merged.top.phylum) # 11 taxa remain, as expected
taxa_names(merged.top.phylum)

# Extract phylum counts & sample metadata into a data frame
phydata<-as(otu_table(merged.top.phylum),'matrix') %>% t() %>% as.data.frame() %>% 
  mutate("sampleid_2"=factor(row.names(.))) %>% 
  merge(.,as(sample_data(merged.top.phylum),'data.frame')%>%  dplyr::select(sampleid_2,OriginalSampleID, TLR3, BirdID, AgeClass5,SexEstimate),by='sampleid_2') %>% 
  gather(key=Phylum,value=Abundance,2:12) %>% 
  group_by(OriginalSampleID) %>% mutate(RelAbund=Abundance/sum(Abundance)) %>%
  ungroup %>% as.data.frame %>%
  mutate(Phylum=ifelse(Phylum=="2480e81716ba1920a54162ee3edbcc99",'Low abundance',Phylum))%>% mutate(Phylum=ordered(Phylum,levels=rev(c('Low abundance',top10phyla[10:1]))))  
view(phydata)

#to get phylum names
write.table(phydata, "phydata.txt", sep="\t", col.names=TRUE, row.names=FALSE)
View(merged.top.phylum@tax_table@.Data)

# colour scheme
phyPalette<-c("dodgerblue",
              "firebrick1",
              "gold",
              "darkorchid",
              "chartreuse4",
              "aquamarine",
              "darkorange1",
              "darkslategrey",
              "deeppink1",
              "greenyellow",
              "grey80")

#import dataframe (with additional order column) back into R - see example file
phydata <- read_csv("phylumtaxaplot.csv", col_types = cols(Order = col_number(),
                                                           Abundance = col_number(),
                                                           RelAbund = col_number()))
View(phydata)


#make phylum graph
taxaphylum<-phydata %>% mutate(AgeClass5=ordered(AgeClass5,levels=(c("CH","FL","OFL","SA","A")))) %>% arrange(AgeClass5) %>%
  mutate(Phylum=ordered(Phylum,levels=(c("Proteobacteria (42%)","Firmicutes (22%)","Actinobacteria (17%)","Chloroflexi (5%)","Planctomycetes (3%)","Bacteroidetes (3%)", "Patescibacteria (2%)","Acidobacteria (1%)","Verrucomicrobia (1%)","Cyanobacteria (1%)","Low abundance (1%)")))) %>% arrange(Phylum) %>%
  ggplot(.,aes(x=Order,y=RelAbund*100,fill=Phylum,color=Phylum)) +
  geom_bar(colour = alpha("white", 0.4), width=1,lwd=0.001, stat='identity') +
  facet_grid(~ AgeClass5, scales = "free", space ="free") +
  theme_classic()+
  ylab("Relative abundance (%)")+
  xlab("Sample")+
  scale_fill_manual(values=phyPalette[1:12])+ guides(fill=guide_legend(reverse=FALSE))+
  theme(legend.title=element_text(size=15,face='bold'),legend.text=element_text(size=10,face='bold'))+  scale_color_manual(values=phyPalette[1:12],guide=FALSE)+
  scale_x_discrete() + scale_y_continuous(expand = c(0, 0))+
  theme(plot.title=element_text(size=10,face='bold'))+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=15,face='bold',vjust=1.5))+
  theme(axis.title.x=element_text(size=15,face='bold'))+
  theme(axis.title.y=element_text(size=15,face='bold'),axis.text.y=element_text(size =10,face='bold'))+ theme(panel.grid.major= element_blank())
plot(taxaphylum)

#B:class graph####

#this gets down into merged columns to make tidier graph have additionally re-ordered in excel

class_merge <- tax_glom(Cousin_final, taxrank="Class") 
View(class_merge@tax_table@.Data)

# Find top 10 phyla by relative abundance
class_relabund  = transform_sample_counts(class_merge, function(x) x / sum(x) )
top10class <- sort(taxa_sums(class_relabund),decreasing=TRUE)[1:10] %>% names() 
rareclass <- setdiff(taxa_names(class_merge),top10class) 
rareclass
top10class

# Merge rare phyla into "Low abundance" category
merged.top.class<-merge_taxa(class_merge,rareclass,archetype="2480e81716ba1920a54162ee3edbcc99") %>%
  prune_taxa(taxa_sums(.)>0,.)
ntaxa(merged.top.class)
taxa_names(merged.top.class)

# Extract class counts & sample metadata into a data frame
classdata<-as(otu_table(merged.top.class),'matrix') %>% t() %>% as.data.frame() %>% 
  mutate("sampleid_2"=factor(row.names(.))) %>% 
  merge(.,as(sample_data(merged.top.class),'data.frame')%>%select(sampleid_2,OriginalSampleID, TLR3, BirdID, AgeClass5,SexEstimate),by='sampleid_2') %>% 
  gather(key=Class,value=Abundance,2:12) %>% 
  group_by(OriginalSampleID) %>% mutate(RelAbund=Abundance/sum(Abundance)) %>%
  ungroup %>% as.data.frame %>%
  mutate(Class=ifelse(Class=="2480e81716ba1920a54162ee3edbcc99",'Low abundance',Class))%>% mutate(Class=ordered(Class,levels=rev(c('Low abundance',top10class[10:1]))))  
view(classdata)

#to get class names
write.table(classdata, "classdata.txt", sep="\t", col.names=TRUE, row.names=FALSE)
View(merged.top.class@tax_table@.Data)

# colour scheme
classPalette<-c("dodgerblue",
                "#a5d2ff",
                "gold",
              "firebrick1",
              "#ffc0c0",
              "darkorchid",
              "chartreuse4",
              "aquamarine",
                "orange",
              "deeppink1",
              "grey80")

#import back into R (with additional order column)
classdata <- read_csv("classtaxaplot.csv", col_types = cols(Abundance = col_number(),
                                                           RelAbund = col_number()))
View(classdata)

#make graph
taxaclass<-classdata %>% mutate(AgeClass5=ordered(AgeClass5,levels=(c("CH","FL","OFL","SA","A")))) %>% arrange(AgeClass5) %>%
  mutate(Class=ordered(Class,levels=(c("Gammaproteobacteria (25%)","Alphaproteobacteria (16%)","Actinobacteria (16%)","Bacilli (16%)","Clostridia (6%)","Chloroflexia (4%)", "Planctomycetacia (3%)","Bacteroidia (3%)","Saccharimonadia (2%)","Verrucomicrobiae (1%)","Low abundance (8%)")))) %>% arrange(Class) %>%
  ggplot(.,aes(x=order,y=RelAbund*100,fill=Class,color=Class)) +
  geom_bar(colour = alpha("white", 0.4), width=1,lwd=0.001, stat='identity') +
  facet_grid(~ AgeClass5, scales = "free", space ="free") +
  theme_classic()+
  ylab("Relative abundance (%)")+
  xlab("Sample")+
  scale_fill_manual(values=classPalette[1:12])+ guides(fill=guide_legend(reverse=FALSE))+
  theme(legend.title=element_text(size=15,face='bold'),legend.text=element_text(size=10,face='bold'))+  scale_color_manual(values=classPalette[1:12],guide=FALSE)+
  scale_x_discrete() + scale_y_continuous(expand = c(0, 0))+
  theme(plot.title=element_text(size=10,face='bold'))+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=15,face='bold',vjust=1.5))+
  theme(axis.title.x=element_text(size=15,face='bold'))+
  theme(axis.title.y=element_text(size=15,face='bold'),axis.text.y=element_text(size =10,face='bold'))+ theme(panel.grid.major= element_blank())
plot(taxaclass)

# add in labels
phylaplot <- arrangeGrob(taxaphylum, top = textGrob("A",  x = unit(0, "npc"), y   = unit(0.4, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=10)))
classplot <- arrangeGrob(taxaclass, top = textGrob("B",  x = unit(0, "npc"), y   = unit(0.4, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=10)))

#export graph
combo_plot<-grid.arrange(phylaplot,classplot, ncol = 1)
ggsave("Fig1.jpeg", combo_plot, width = 18, height =12 ,units="cm", dpi = 300)



##################################### Extract core Microbiome ###########################

# using microbiome####
library(microbiome)

#to plot heatmap of core microbiome - family

library(knitr)
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
#

family_merge <- tax_glom(Cousin_final, taxrank="Family")
#relative abundance
pseq.rel <- microbiome::transform(family_merge, "compositional")

# Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(pseq.rel, detection = 0.1/100, sort = TRUE))

#plot core line plot
plot_core(pseq.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")


familycore <- plot_core(pseq.rel, plot.type = "heatmap",
                        prevalences = prevalences,
                        detections = detections,
                        min.prevalence = .1, horizontal = F)+
  xlab("Detection Threshold (Relative Abundance (%))")

# get the data used for plotting 
df <- familycore$data
# get the list of OTUs
list <- df$Taxa
# check the OTU ids
# print(list)
#this tells you the prevelance of each merged taxa - but it does so for various different detection thresholds - hence why there are ~10 outputs for each taxa
View(df)

#this gives name of top ASVs
core.taxa.standard <- core_members(pseq.rel, detection = 0.001, prevalence = 0.5)

# gives you the taxa name for taxa in at least 50% of samples

taxonomy <- as.data.frame(tax_table(pseq.rel))

core.taxa.standard <- core_members(pseq.rel, detection = 0.001, prevalence = 50/100)

#print(core.taxa.standard)

# Subset this taxonomy table to include only core OTUs  
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard)

DT::datatable(core_taxa_id)

# create table S4 using df and core_taxa_id


# repeat but using genus instead of family ####

genus_merge <- tax_glom(Cousin_final, taxrank="Genus")

#relative abundance
pseq.rel.gen <- microbiome::transform(genus_merge, "compositional")

# Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(pseq.rel.gen, detection = 0.1/100, sort = TRUE))

#plot core line plot
plot_core(pseq.rel.gen, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")


genuscore <- plot_core(pseq.rel.gen, plot.type = "heatmap",
                       prevalences = prevalences,
                       detections = detections,
                       min.prevalence = .1, horizontal = F)+
  xlab("Detection Threshold (Relative Abundance (%))")

# get the data used for plotting 
df <- genuscore$data
# get the list of OTUs
list <- df$Taxa
# check the OTU ids
# print(list)

#this tells you the prevelance of each merged taxa - but it does so for various different detection thresholds - hence why there are ~10 outputs for each taxa
View(df)

#this gives name of top ASVs
core.taxa.standard.gen <- core_members(pseq.rel.gen, detection = 0.001, prevalence = 0.5)

# gives you the taxa name for taxa in at least 80% of samples

taxonomy <- as.data.frame(tax_table(pseq.rel.gen))

core.taxa.standard.gen <- core_members(pseq.rel.gen, detection = 0.001, prevalence = 50/100)

#print(core.taxa.standard)

# Subset this taxonomy table to include only core OTUs  
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard.gen)

DT::datatable(core_taxa_id)


##################################### GM repeatability -  Fig S3 and S4 #########
# Individual repeatability ####

#For duplicated birds (not including duplicated sequencing)

# Subset dataset

# removing duplicated samples (sequencing repeats)
Brepeat <- subset_samples(Cousin_final , Keep%in%c("y"))
# Keep only birds repeated within a season
BirdRepeat2 <- subset_samples(Brepeat, BirdSeasonRepeat%in%c("Yes"))
# only ASVs with read depth greater than 50
BirdRepeat <-prune_taxa(taxa_sums(BirdRepeat2 )>50,BirdRepeat2)

# get dissimilarity matrix for repeat samples 

# extract alpha diversity metrics of choice from phyloseq then import into a table

#alpha diversity

set.seed(28367)
rare_birdrepeat = rarefy_even_depth(BirdRepeat, sample.size = 10000)
birdrepeat_rare.alpha.diversity <- estimate_richness(rare_birdrepeat, measures = c( "Shannon"))

# calculate faiths PD
birdrepeat_rare_faith_PD <- estimate_pd(rare_birdrepeat)

# bind tigether
birdrepeat_rarealphadata <- cbind(sample_data(rare_birdrepeat), birdrepeat_rare.alpha.diversity, birdrepeat_rare_faith_PD)

write.table(birdrepeat_rarealphadata, "birdrepeat_rarealphadata.txt")

# extract beta dissimilarity
library(vegan)

## Further filtering: Keep taxa when appearing in minimum % samples
prev0 <- apply(X = otu_table(BirdRepeat),
               MARGIN = ifelse(taxa_are_rows(BirdRepeat), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# taxa in 2.5% of samples
BirdRepeat_clean <- prune_taxa((prev0 > 5), BirdRepeat)
ntaxa(BirdRepeat_clean) # 1011

weighted <- diversity_comparison(BirdRepeat_clean, distance = "wunifrac")
unweighted <- diversity_comparison(BirdRepeat_clean, distance = "unifrac")

#export tables for each diversity metric of interest 
write.table(weighted, "birdweighted.txt", sep="\t", col.names=TRUE, row.names=FALSE)
write.table(unweighted, "birdunweighted.txt", sep="\t", col.names=TRUE, row.names=FALSE)
write.table(shannon, "birdshannon.txt", sep="\t", col.names=TRUE, row.names=FALSE)

# In excel seperate within vs between samples, remove duplicated values and calculate euclidean distances

#import back into R
bird_repeated <- read_csv("repeatedbird_rare_distance.csv")
bird_repeated$Weighted <-as.numeric(bird_repeated$weighted)
bird_repeated$UnWeighted <-as.numeric(bird_repeated$unweighted)
bird_repeated$shannon_difference <-as.numeric(bird_repeated$Shannon_difference)
view(bird_repeated)

#remove the two shannon diversity outliers
bird_repeated_noout<- bird_repeated %>%  filter(Remove_outlier != "Yes")

# Run models

shannon_catch_rarealpha_birdrepeated_k<-kruskal.test(shannon_difference~BirdRepeat, data=bird_repeated)
shannon_catch_rarealpha_birdrepeated_k
dunnTest(shannon_difference~CatchRepeat, method="bh",data=bird_repeated)
shannon_catch_rarealpha_birdrepeated_k<-kruskal.test(shannon_difference~BirdRepeat, data=bird_repeated_noout)
shannon_catch_rarealpha_birdrepeated_k
dunnTest(shannon_difference~CatchRepeat, method="bh",data=bird_repeated_noout)

w_catch_rarealpha_birdrepeated_k<-kruskal.test(Weighted~BirdRepeat, data=bird_repeated)
w_catch_rarealpha_birdrepeated_k
dunnTest(Weighted~CatchRepeat, method="bh",data=bird_repeated)

unw_catch_rarealpha_birdrepeated_k<-kruskal.test(UnWeighted~BirdRepeat, data=bird_repeated)
unw_catch_rarealpha_birdrepeated_k
dunnTest(UnWeighted~CatchRepeat,method="bh", data=bird_repeated)

# Fig S3

BR_plot1<-ggplot(bird_repeated, aes(x = BirdRepeat, y = shannon_difference)) +theme_classic()+geom_jitter(aes(colour = Remove_outlier),alpha=0.5, position = position_jitter(0.2))+geom_boxplot(outlier.shape = NA, size = .75, width=0.75, alpha = 0.4)+stat_summary(fun = "mean", geom = "point", shape = 5, size = 2, color = "grey32")+labs( y="Shannon dissimilarity") + theme(axis.title.x = element_blank(),axis.text.x =element_text(size=18))+  theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=18)) +scale_x_discrete(labels=c("Different individual","Same individual"))+theme(legend.position = "none")+scale_color_manual( values = c("lightgrey","#0B5F6F" ))

BR_plot2<-ggplot(bird_repeated, aes(x = BirdRepeat, y = UnWeighted)) +theme_classic()+geom_jitter(alpha=0.5, position = position_jitter(0.2), color="lightgrey")+geom_boxplot(outlier.shape = NA, size = .75, width=0.75, alpha = 0.4)+stat_summary(fun = "mean", geom = "point", shape = 5, size = 2, color = "grey32")+labs( y="Unweighted unifrac disimilarity") + theme(axis.title.x = element_blank(),axis.text.x =element_text(size=18))+  theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=18)) +scale_x_discrete(labels=c("Different individual","Same individual"))

BR_plot3<-ggplot(bird_repeated, aes(x = BirdRepeat, y = Weighted)) +theme_classic()+geom_jitter(alpha=0.5, position = position_jitter(0.2), color="lightgrey")+geom_boxplot(outlier.shape = NA, size = .75, width=0.75, alpha = 0.4)+stat_summary(fun = "mean", geom = "point", shape = 5, size = 2, color = "grey32")+labs( y="Weighted unifrac disimilarity") +  theme(axis.title.x = element_blank(),axis.text.x =element_text(size=18))+  theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=18)) +scale_x_discrete(labels=c("Different individual","Same individual"))

library(gridExtra)

BR_plot1 <- arrangeGrob(BR_plot1, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

BR_plot2 <- arrangeGrob(BR_plot2, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

BR_plot3 <- arrangeGrob(BR_plot3, top = textGrob("C", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

combo_plot<-grid.arrange( BR_plot1,BR_plot2,BR_plot3, ncol = 3)
ggsave("Fig_S3.png", combo_plot, width = 18, height = 7, dpi = 300)

# Sequencing repeatability####
# using duplicated sequencing samples - previously filtered and saved as Repeat phyloseq object

# Fig S4A - make taxa graph of repeated samples only
# this gets down into merged collumns but need to re-order slighty in excell.....

repeat_phylum_merge <- tax_glom(Repeat, taxrank="Phylum") 
View(repeat_phylum_merge@tax_table@.Data)

# Find top 10 phyla by abundance
sort(taxa_sums(repeat_phylum_merge),decreasing=TRUE)[1:10] %>% names() -> top10phyla
setdiff(taxa_names(repeat_phylum_merge),top10phyla) -> rarephyla
rarephyla
top10phyla

# Merge rare phyla into "Low abundance" category
merged.top.phylum<-merge_taxa(repeat_phylum_merge,rarephyla,archetype="f31d538f18b62ff69bc47fc23b7c4f81") %>%
  prune_taxa(taxa_sums(.)>0,.)
ntaxa(merged.top.phylum)
taxa_names(merged.top.phylum)

# Extract phylum counts & sample metadata into a data frame
repeat_phydata<-as(otu_table(merged.top.phylum),'matrix') %>% t() %>% as.data.frame() %>% 
  mutate("sampleid_2"=factor(row.names(.))) %>% 
  merge(.,as(sample_data(merged.top.phylum),'data.frame')%>%select(sampleid_2,OrigiNAlSampleID,OriginalSampleIDUnique, TLR3, BirdID, AgeClass4,SexEstimate, CleanedFeatureCount, Keep_Catch_Bag_Tray),by='sampleid_2') %>% 
  gather(key=Phylum,value=Abundance,2:12) %>% 
  group_by(OriginalSampleIDUnique) %>% mutate(RelAbund=Abundance/sum(Abundance)) %>%  
  ungroup %>% as.data.frame %>%
  mutate(Phylum=ifelse(Phylum=="f31d538f18b62ff69bc47fc23b7c4f81",'Low abundance',Phylum))%>% mutate(Phylum=ordered(Phylum,levels=rev(c('Low abundance',top10phyla[10:1]))))  
view(repeat_phydata)

# to get class names
write.table(repeat_phydata, "repeat_phydata.txt", sep="\t", col.names=TRUE, row.names=FALSE)
View(merged.top.phylum@tax_table@.Data)

# colour scheme
phyPalette<-c("dodgerblue",
              "firebrick1",
              "gold",
              "darkorchid",
              "chartreuse4",
              "orange",
              "aquamarine",
              "darkslategrey",
              "deeppink1",
              "greenyellow",
              "grey80")

# import back into R (with additional order column)
repeat_phydata <- read_csv("repeat_phylumtaxaplot.csv")
View(repeat_phydata)


# make graph  - Fig S4A
repeat_taxaphylum<-repeat_phydata %>% mutate(AgeClass4=ordered(AgeClass4,levels=(c("CH","FL","SA","A")))) %>% arrange(AgeClass4) %>%
  mutate(Phylum=ordered(Phylum,levels=(c("Proteobacteria","Firmicutes","Actinobacteria","Chloroflexi","Planctomycetes","Bacteroidetes", "Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Low abundance")))) %>% arrange(Phylum) %>%
  ggplot(.,aes(x=Order,y=RelAbund*100,fill=Phylum,color=Phylum)) +
  geom_bar(colour = alpha("white", 0.4), width=1,lwd=0.001, stat='identity') +
  theme_classic()+
  ylab("Relative abundance (%)")+
  xlab("Sample")+
  scale_fill_manual(values=phyPalette[1:12])+ guides(fill=guide_legend(reverse=FALSE))+
  theme(legend.title=element_text(size=15,face='bold'),legend.text=element_text(size=10,face='bold'))+  scale_color_manual(values=phyPalette[1:12],guide=FALSE)+
  scale_x_discrete() + scale_y_continuous(expand = c(0, 0))+
  theme(plot.title=element_text(size=10,face='bold'))+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=15,face='bold',vjust=1.5))+
  theme(axis.title.x=element_text(size=15,face='bold'))+
  theme(axis.title.y=element_text(size=15,face='bold'),axis.text.y=element_text(size =10,face='bold'))+ theme(panel.grid.major 
                                                                                                              = element_blank())
ggsave("Fig_S4a.png", repeat_taxaphylum, width = 18, height = 7, dpi = 300)

# pairwise Euclidean dissimilarity between different samples, versus within pairs of duplicated samples 

# extract alpha diversity metrics of choice from phyloseq then import into a table

#repeat and alpha diversity

set.seed(28367)
rare_repeat = rarefy_even_depth(Repeat, sample.size = 10000)
repeat_rarealphadata <- estimate_richness(rare_repeat, measures = c( "Shannon"))

write.table(repeat_rarealphadata, "repeat_rarealphadata.txt")

#read table back in
repeated_alpha <- read_csv("repeated_alpha_distance.csv", col_types = cols(Shannon_difference = col_number()))

#make graph Fig S4B
SR_plot1<- ggplot(repeated_alpha, aes(x = between_within, y = Shannon_difference)) +theme_classic()+geom_jitter(alpha=0.5, position = position_jitter(0.2), color="lightgrey")+geom_boxplot(outlier.shape = NA, size = .75, width=0.75, alpha = 0.4)+stat_summary(fun = "mean", geom = "point", shape = 5, size = 2, color = "grey32")+labs(x="Sample", y="Rrefied shannon dissimilarity)") + theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=18))+  theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=18)) +scale_x_discrete(labels=c("between","within"))

# beta dissimilarity
library(vegan)

weighted <- diversity_comparison(Repeat, distance = "wunifrac")
unweighted <- diversity_comparison(Repeat, distance = "unifrac")

#export tables and seperate within vs between
write.table(weighted, "weighted.txt", sep="\t", col.names=TRUE, row.names=FALSE)
write.table(unweighted, "unweighted.txt", sep="\t", col.names=TRUE, row.names=FALSE)

#import back into R
repeated <- read_csv("REPEATED_MEASURES_COMBINED.csv", col_types = cols(unweighted = col_number(),weighted = col_number()))
view(repeated)

# remove samples that no longer have a duplicate
repeated_measures  <- repeated %>% filter(keep_remove != "REMOVE")

SR_plot2<ggplot(repeated_measures, aes(x = between_within, y = unweighted)) +theme_classic()+geom_jitter(alpha=0.5, position = position_jitter(0.2), color="lightgrey")+geom_boxplot(outlier.shape = NA, size = .75, width=0.75, alpha = 0.4)+stat_summary(fun = "mean", geom = "point", shape = 5, size = 2, color = "grey32")+labs(x="Sample", y="Unweighted unifrac disimilarity)") + theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=18))+  theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=18)) +scale_x_discrete(labels=c("between","within"))

SR_plot3<ggplot(repeated_measures, aes(x = between_within, y = weighted)) +theme_classic()+geom_jitter(alpha=0.5, position = position_jitter(0.2), color="lightgrey")+geom_boxplot(outlier.shape = NA, size = .75, width=0.75, alpha = 0.4)+stat_summary(fun = "mean", geom = "point", shape = 5, size = 2, color = "grey32")+labs(x="Sample", y="Weighted unifrac disimilarity)") + theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=18))+  theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=18)) +scale_x_discrete(labels=c("between","within"))

# Run models

shannon_catch_rarealpha_repeated_k<-kruskal.test(shannon_difference~between_within, data=repeated_alpha)
shannon_catch_rarealpha_repeated_k
dunnTest(shannon_difference~between_within, method="bh",data=repeated_alpha)

w_catch_rarealpha_repeated_k<-kruskal.test(Weighted~between_within, data=repeated_measures)
w_catch_rarealpha_repeated_k
dunnTest(Weighted~between_within, method="bh",data=repeated_measures)

unw_catch_rarealpha_repeated_k<-kruskal.test(unweighted~between_within, data=repeated_measures)
unw_catch_rarealpha_repeated_k
dunnTest(unweighted~between_within, method="bh",data=repeated_measures)

library(gridExtra)

SR_plot1 <- arrangeGrob(SR_plot1, top = textGrob("i", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

SR_plot2 <- arrangeGrob(SR_plot2, top = textGrob("ii", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

SR_plot3 <- arrangeGrob(SR_plot3, top = textGrob("iii", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

combo_plot<-grid.arrange( SR_plot1,SR_plot2,SR_plot3, ncol = 3)
ggsave("Fig_S4.png", combo_plot, width = 18, height = 7, dpi = 300)# Add Fig S4a to plot

#####################################ALPA DIVERSITY ANALYSIS #############################################


# use Cousin_final_cleaned and Cousin_final_cleaned_rare DATASETS

#extract alpha diversity metrics of choice from phyloseq then import into a table####
alpha.diversity <- estimate_richness(Cousin_final_cleaned, measures = c( "Chao1", "Shannon"))
rare.alpha.diversity <- estimate_richness(Cousin_final_cleaned_rare, measures = c( "Chao1", "Shannon"))

# calculate faiths PD
faith_PD <- estimate_pd(Cousin_final_cleaned)
rare_faith_PD <- estimate_pd(Cousin_final_cleaned_rare)

# bind tigether
head(alpha.diversity)
alphadata <- cbind(sample_data(Cousin_final_cleaned), alpha.diversity, faith_PD)
rarealphadata <- cbind(sample_data(Cousin_final_cleaned_rare), rare.alpha.diversity, rare_faith_PD)
View(alphadata)
View(rarealphadata)

#write.table(alphadata, "alphadiversity.txt")
#write.table(rarealphadata, "rarealphadiversity.txt")

#remove samples with no MHC  info - alpha data
alphadata_all<- alphadata %>%  filter(Both_MHC_present != "N")

#make sure variables coded properly
alphadata_all$BirdID <- as.factor(alphadata_all$BirdID)
alphadata_all$RunID <- as.factor(alphadata_all$RunID)
alphadata_all$CatchID <- as.factor(alphadata_all$CatchID)
alphadata_all$Ase.dab3 <- as.factor(alphadata_all$Ase.dab3)
alphadata_all$Ase.dab4 <- as.factor(alphadata_all$Ase.dab4)
alphadata_all$Ase.dab5 <- as.factor(alphadata_all$Ase.dab5)
alphadata_all$Ase-ua1 <- as.factor(alphadata_all$Ase-ua1)
alphadata_all$Ase-ua3 <- as.factor(alphadata_all$Ase-ua3)
alphadata_all$Ase-ua4 <- as.factor(alphadata_all$Ase-ua4)
alphadata_all$Ase-ua5 <- as.factor(alphadata_all$Ase-ua5)
alphadata_all$Ase-ua7 <- as.factor(alphadata_all$Ase-ua7)
alphadata_all$Ase-ua8 <- as.factor(alphadata_all$Ase-ua8)
alphadata_all$Ase-ua6 <- as.factor(alphadata_all$Ase-ua6)
alphadata_all$Ase-ua9 <- as.factor(alphadata_all$Ase-ua9)
alphadata_all$MAse-ua11 <- as.factor(alphadata_all$Ase-ua11)
alphadata_all$MAse-ua10 <- as.factor(alphadata_all$Ase-ua10)
alphadata_all$RunID <- as.factor(alphadata_all$RunID)

#remove samples with no MHC  info - rarefied alpha data and check factors are factors
rarealphadata_all<- rarealphadata %>%  filter(Both_MHC_present != "N")

#make sure variables coded properly
rarealphadata_all$BirdID <- as.factor(rarealphadata_all$BirdID)
rarealphadata_all$RunID <- as.factor(rarealphadata_all$RunID)
rarealphadata_all$CatchID <- as.factor(rarealphadata_all$CatchID)
rarealphadata_all$Ase-dab3 <- as.factor(rarealphadata_all$Ase-dab3)
rarealphadata_all$Ase-dab4 <- as.factor(rarealphadata_all$Ase-dab4)
rarealphadata_all$Ase-dab5 <- as.factor(rarealphadata_all$Ase-dab5)
rarealphadata_all$Ase-ua1 <- as.factor(rarealphadata_all$Ase-ua1)
rarealphadata_all$Ase-ua3 <- as.factor(rarealphadata_all$Ase-ua3)
rarealphadata_all$Ase-ua4 <- as.factor(rarealphadata_all$Ase-ua4)
rarealphadata_all$Ase-ua5 <- as.factor(rarealphadata_all$Ase-ua5)
rarealphadata_all$Ase-ua7 <- as.factor(rarealphadata_all$Ase-ua7)
rarealphadata_all$Ase-ua8 <- as.factor(rarealphadata_all$Ase-ua8)
rarealphadata_all$Ase-ua6 <- as.factor(rarealphadata_all$Ase-ua6)
rarealphadata_all$Ase-ua9 <- as.factor(rarealphadata_all$Ase-ua9)
rarealphadata_all$Ase-ua11 <- as.factor(rarealphadata_all$Ase-ua11)
rarealphadata_all$Ase-ua10 <- as.factor(rarealphadata_all$Ase-ua10)
rarealphadata_all$RunID <- as.factor(rarealphadata_all$RunID)

#remove chick category
alphadata_allc<- alphadata_all %>% filter(AgeClass5 != "CH")
rarealphadata_allc<- rarealphadata_all %>% filter(AgeClass5 != "CH")

#remove duplicated birds  - chosen at random 
alphadata_allc_nodup<- alphadata_allc %>% filter(Keep_DUPLICATE_BIRD != "n")
rarealphadata_allc_nodup<- rarealphadata_allc %>% filter(Keep_DUPLICATE_BIRD != "n")
view(rarealphadata_allc_nodup)

#check normality - all are not normal - check model outputs to determine next steps etc- ####
shapiro.test(alphadata_allc_nodup$Shannon)
shapiro.test(alphadata_allc_nodup$Chao1)
shapiro.test(alphadata_allc_nodup$PD)
shapiro.test(rarealphadata_allc_nodup$Shannon)
shapiro.test(rarealphadata_allc_nodup$Chao1)
shapiro.test(rarealphadata_allc_nodup$PD)

# scale and centralise variables ####

#load packages
library(MuMIn)
library(DHARMa)
library(sjPlot)
library(arm)
library(visreg)
library(car)

#scale and centralise continuos varaibles for chick and duplicate removed dataset 
rarealphadata_allc_nodup$ScaledHs_obs<-rescale(rarealphadata_allc_nodup$Hs_obs, binary.inputs="full")
alphadata_allc_nodup$ScaledHs_obs<-rescale(alphadata_allc_nodup$Hs_obs, binary.inputs="full")

rarealphadata_allc_nodup$ScaledMHC1_Diversity<-rescale(rarealphadata_allc_nodup$MHC1_Diversity, binary.inputs="full")
alphadata_allc_nodup$ScaledMHC1_Diversity<-rescale(alphadata_allc_nodup$MHC1_Diversity, binary.inputs="full")
rarealphadata_allc_nodup$ScaledMHC1_Diversity_2<-rescale(rarealphadata_allc_nodup$MHC1_Diversity_2, binary.inputs="full")
alphadata_allc_nodup$ScaledMHC1_Diversity_2<-rescale(alphadata_allc_nodup$MHC1_Diversity_2, binary.inputs="full")

rarealphadata_allc_nodup$ScaledMHC2_Diversity<-rescale(rarealphadata_allc_nodup$MHC2_Diversity, binary.inputs="full")
alphadata_allc_nodup$ScaledMHC2_Diversity<-rescale(alphadata_allc_nodup$MHC2_Diversity, binary.inputs="full")
rarealphadata_allc_nodup$ScaledMHC2_Diversity_2<-rescale(rarealphadata_allc_nodup$MHC2_Diversity_2, binary.inputs="full")
alphadata_allc_nodup$ScaledMHC2_Diversity_2<-rescale(alphadata_allc_nodup$MHC2_Diversity_2, binary.inputs="full")

# explore rarefied dataset ####
par(mfrow=c(2,2))
hist(rarealphadata_allc_nodup$Shannon)
hist(log(rarealphadata_allc_nodup$Shannon+1))
hist(log10(rarealphadata_allc_nodup$Shannon))
hist(sqrt(rarealphadata_allc_nodup$Shannon))
par(mfrow=c(1,1)) 
#shannon looks okay not transformed - little skewed
par(mfrow=c(2,2))
hist(rarealphadata_allc_nodup$Chao1)
hist(log(rarealphadata_allc_nodup$Chao1+1))
hist(log10(rarealphadata_allc_nodup$Chao1))
hist(sqrt(rarealphadata_allc_nodup$Chao1))
par(mfrow=c(1,1)) 
#chao1 may be better log or sqrt transformed - see model residuals
par(mfrow=c(2,2))
hist(rarealphadata_allc_nodup$PD)
hist(log(rarealphadata_allc_nodup$PD+1))
hist(log10(rarealphadata_allc_nodup$PD))
hist(sqrt(rarealphadata_allc_nodup$PD))
par(mfrow=c(1,1)) 
#PD may be better log or sqrt transformed - see model residuals

#make sure no difference between runs
shannon_run<-lm(Shannon~RunID, data=rarealphadata_allc_nodup)
Anova(shannon_run)
chao_run<-lm(Chao1~RunID, data=rarealphadata_allc_nodup)
Anova(chao_run)
faith_run<-lm(PD~RunID, data=rarealphadata_allc_nodup)
Anova(faith_run)

#make sure no difference between AM/PM
shannon_time<-lm(Shannon~CollectionTime, data=alphadata_all)
Anova(shannon_time)
chao_time<-lm(Chao1~CollectionTime, data=alphadata_all)
Anova(chao_time)
faith_time<-lm(PD~CollectionTime, data=alphadata_all)
Anova(faith_time)

# check for overall collineairty between variables
library("psych")
scatterplotMatrix(~AgeClass5+AgeClass5+SexEstimate+TLR3+ScaledHs_obs+ScaledMHC1_Diversity+ScaledMHC2_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity_2+RunID+FieldPeriodID+TQ, data=rarealphadata_allc_nodup)
# runID and field period a little correlated

#test whether presence/absence of MHC alleles are correlated
library(GGally)

MHC_allele_CORRELATION<-ggpairs(rarealphadata_allc_nodup,columns=c("TLR3","Ase-ua5","Ase-ua3","Ase-ua6","Ase-ua7","Ase-ua9","Ase-ua11","Ase-ua8","Ase-ua4", "Ase-ua1","Ase-ua10", "Ase-dab3","Ase-dab4","Ase-dab5"))
MHC_allele_CORRELATION

# convert factors to numeric
cor_mhc<- rarealphadata_allc_nodup[ , c("TLR3","Ase-ua5","Ase-ua3","Ase-ua6","Ase-ua7","Ase-ua9","Ase-ua11","Ase-ua8","Ase-ua4", "Ase-ua1", "Ase-ua10","Ase-dab3","Ase-dab4","Ase-dab5")]
MHC_allele_CORREL<-ggcorr(cor_mhc,label=TRUE)
MHC_allele_CORREL


################## First test for differences in diversity - using Cousin_final_cleaned or Cousin_final_cleaned_rare ####

#Shannon full model - non-rarefied - Diversity####


#full model 
shannon2_nodup <- lm(Shannon ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2+CleanedFeatureCount, data = alphadata_allc_nodup , na.action = "na.fail")

#plot residuaals and check effects 
resp = simulateResiduals(shannon2_nodup )
plot(resp, rank = T)
plot_model(shannon2_nodup )

# to test if residuals are underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = shannon2_nodup , n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(shannon2_nodup )
Anova(shannon2_nodup )
#visreg(shannon2_nodup)
model_shannon2_nodup <-plot_model(shannon2_nodup )
summ(shannon2_nodup )

#use MUMIN::DREDGE to test all models options and find best one and average
shannon2_nodup_dredge<-dredge(shannon2_nodup)
shannon2_nodup_dredge
model.avg(shannon2_nodup_dredge)
#plot(shannon2_dredge)
averaged_shannon2_nodup<-summary(model.avg(get.models(shannon2_nodup_dredge, subset = delta < 7)))
averaged_shannon2_nodup
sw(averaged_shannon2_nodup)

#chao1 full model - non-rarefied - Diversity####
chao2a_nodup <- lm(Chao1 ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

#plots look bad - try log transforming
alphadata_allc_nodup$logchao <- log((alphadata_allc_nodup$Chao1)+1)

chao2_nodup <- lm(logchao ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

#plot residuaals and check effects - lookls fine when log transformed
resp = simulateResiduals(chao2a_nodup)
plot(resp, rank = T)
resp = simulateResiduals(chao2_nodup)
plot(resp, rank = T)
#visreg(rare_chao2_nodup)
plot_model(chao2_nodup)
model_chao2_nodup<-plot_model(chao2_nodup)

# to test if residuals are underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = chao2_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1))

#use MUMIN::DREDGE to test all models options and find best one and average
chao2_nodup_dredge<-dredge(chao2_nodup)
chao2_nodup_dredge
model.avg(chao2_nodup_dredge)
#plot(chao2_dredge)
averaged_chao2_nodup<-summary(model.avg(get.models(chao2_nodup_dredge, subset = delta < 7)))
averaged_chao2_nodup
sw(averaged_chao2_nodup)

#faith full model - non-rarefied - diversity ####
Faith2a_nodup <- lm(PD ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

resp = simulateResiduals(Faith2a_nodup)
plot(resp, rank = T) 

#plots look bad - try log transforming
alphadata_allc_nodup$logFaith <- log((alphadata_allc_nodup$PD)+1)

Faith2_nodup <- lm(logFaith ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

#plot residuals and check effects - look alot better using log
resp = simulateResiduals(Faith2_nodup)
plot(resp, rank = T) 
plot_model(Faith2_nodup)
model_Rare_faith2_nodup<-plot_model(rare_Faith2_nodup)
#visreg(rare_Faith2_nodup)

# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = Faith2_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1))


#use MUMIN::DREDGE to test all models options and find best one and average
Faith2_nodup_dredge<-dredge(Faith2_nodup)
Faith2_nodup_dredge
model.avg(Faith2_nodup_dredge)
#plot(rare_Faith2_dredge)
averaged_Faith2_nodup<-summary(model.avg(get.models(Faith2_nodup_dredge, subset = delta < 7)))
averaged_Faith2_nodup
sw(averaged_Faith2_nodup)


################### Repeat but using rarefied data ####
#Shannon full model - rarefied - Diversity####


#full model for MHC2
rare_shannon2_nodup <- lm(Shannon ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2, data = rarealphadata_allc_nodup , na.action = "na.fail")

#plot residuals and check effects - look good
resp = simulateResiduals(rare_shannon2_nodup )
plot(resp, rank = T)
plot_model(rare_shannon2_nodup )


# to test if residuals are underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = rare_shannon2_nodup , n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(rare_shannon2_nodup )
Anova(rare_shannon2_nodup )
#visreg(rare_shannon2_nodup)
model_rare_shannon2_nodup <-plot_model(rare_shannon2_nodup )
summ(rare_shannon2_nodup )

#use MUMIN::DREDGE to test all models options and find best one and average
rare_shannon2_nodup_dredge<-dredge(rare_shannon2_nodup)
rare_shannon2_nodup_dredge
model.avg(rare_shannon2_nodup_dredge)
#plot(rare_shannon2_dredge)
averaged_rare_shannon2_nodup<-summary(model.avg(get.models(rare_shannon2_nodup_dredge, subset = delta < 7)))
averaged_rare_shannon2_nodup
sw(averaged_rare_shannon2_nodup)

pairwise.wilcox.test(alphadata_allc$Shannon, alphadata_allc$AgeClass5, p.adjust.method="fdr")

#chao1 full model - rarefied - Diversity####
rare_chao2a_nodup <- lm(Chao1 ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2, data = rarealphadata_allc_nodup, na.action = "na.fail")

#model not a good fit - try log transforming
rarealphadata_allc_nodup$logchao <- log((rarealphadata_allc_nodup$Chao1)+1)

rare_chao2_nodup <- lm(logchao ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2, data = rarealphadata_allc_nodup, na.action = "na.fail")

#plot residuals and check effects - looks fine when log transformed
resp = simulateResiduals(rare_chao2a_nodup)
plot(resp, rank = T)
resp = simulateResiduals(rare_chao2_nodup)
plot(resp, rank = T)
#visreg(rare_chao2_nodup)
plot_model(rare_chao2_nodup)
model_rare_chao2_nodup<-plot_model(rare_chao2_nodup)

# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = rare_chao2_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1))

#use MUMIN::DREDGE to test all models options and find best one and average
rare_chao2_nodup_dredge<-dredge(rare_chao2_nodup)
rare_chao2_nodup_dredge
model.avg(rare_chao2_nodup_dredge)
#plot(rare_chao2_dredge)
averaged_rare_chao2_nodup<-summary(model.avg(get.models(rare_chao2_nodup_dredge, subset = delta < 7)))
averaged_rare_chao2_nodup
sw(averaged_rare_chao2_nodup)

#faith full model - rarefied - Diversity ####
rare_Faith2a_nodup <- lm(PD ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2, data = rarealphadata_allc_nodup, na.action = "na.fail")

resp = simulateResiduals(rare_Faith2a_nodup)
plot(resp, rank = T) 

#model not a good fit - try log transforming
rarealphadata_allc_nodup$logFaith <- log((rarealphadata_allc_nodup$PD)+1)

rare_Faith2_nodup <- lm(logFaith ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+TLR3+ScaledMHC1_Diversity+ScaledMHC1_Diversity_2+ScaledMHC2_Diversity+ScaledMHC2_Diversity_2, data = rarealphadata_allc_nodup, na.action = "na.fail")

#plot residuals and check effects - look alot better using log
resp = simulateResiduals(rare_Faith2_nodup)
plot(resp, rank = T) 
plot_model(rare_Faith2_nodup)
model_Rare_faith2_nodup<-plot_model(rare_Faith2_nodup)
#visreg(rare_Faith2_nodup)

# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = rare_Faith2_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1))


#use MUMIN::DREDGE to test all models options and find best one and average
rare_Faith2_nodup_dredge<-dredge(rare_Faith2_nodup)
rare_Faith2_nodup_dredge
model.avg(rare_Faith2_nodup_dredge)
#plot(rare_Faith2_dredge)
averaged_rare_Faith2_nodup<-summary(model.avg(get.models(rare_Faith2_nodup_dredge, subset = delta < 7)))
averaged_rare_Faith2_nodup
sw(averaged_rare_Faith2_nodup)

######################### Second repeat for Alleles - using non-rarefied data########

#Shannon differences - non-rarefied - alleles####

#looks like no variance in bird ID - remove
shannon1_nodup <- lm(Shannon ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

#check collinearity for each variable
vif(shannon1_nodup)

#plot residuaals and check effects 
resp = simulateResiduals(shannon1_nodup)
plot(resp, rank = T)
plot_model(shannon1_nodup)

# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = shannon1_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(shannon1_nodup)
Anova(shannon1_nodup)
#visreg(rare_shannon1_nodup)
model_shannon1_nodup<-plot_model(shannon1_nodup)
summ(shannon1_nodup)
vif(shannon1_nodup)

#use MUMIN::DREDGE to test all models options and find best one and average
shannon1_nodup_dredge<-dredge(shannon1_nodup)
shannon1_nodup_dredge
model.avg(shannon1_nodup_dredge)
#plot(shannon1_dredge)
averaged_shannon1_nodup<-summary(model.avg(get.models(shannon1_nodup_dredge, subset = delta < 7)))
averaged_shannon1_nodup
sw(averaged_shannon1_nodup)

#chao1 full model - non-rarefied - alleles ####

#full model 
chao1a_nodup <- lm(Chao1 ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5+CleanedFeatureCount, data =alphadata_allc_nodup, na.action = "na.fail")

#model not a good fit - try log transforming
chao1_nodup <- lm(logchao ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

#check collinearity for each variable
vif(chao1_nodup)

#plot residuaals and check effects - look good
resp = simulateResiduals(chao1_nodup)
plot(resp, rank = T)
plot_model(chao1_nodup)
model_chao1_nodup<-plot_model(chao1_nodup)


# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = chao1_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(chao1_nodup)
Anova(chao1_nodup)
#visreg(chao1_nodup)
plot_model(chao1_nodup)
vif(chao1_nodup)

#use MUMIN::DREDGE to test all models options and find best one and average
chao1_nodup_dredge<-dredge(chao1_nodup)
chao1_nodup_dredge
model.avg(chao1_nodup_dredge)
#plot(chao1_dredge)
averaged_chao1_nodup<-summary(model.avg(get.models(chao1_nodup_dredge, subset = delta < 7)))
averaged_chao1_nodup
sw(averaged_chao1_nodup)

#faith full model - non-rarefied - alleles ####

#full model 
faith1a_nodup <- lm(PD ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5+CleanedFeatureCount, data =alphadata_allc_nodup, na.action = "na.fail")

#model not a good fit - try log transforming
faith1_nodup <- lm(logFaith ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5+CleanedFeatureCount, data = alphadata_allc_nodup, na.action = "na.fail")

#check collinearity for each variable
vif(faith1_nodup)

#plot residuals and check effects - look good
resp = simulateResiduals(faith1_nodup)
plot(resp, rank = T)
plot_model(faith1_nodup)
resp = simulateResiduals(faith1a_nodup)
plot(resp, rank = T)
model_faith1<-plot_model(faith1_nodup)


# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = faith1_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(faith1_nodup)
Anova(faith1_nodup)
#visreg(faith1_nodup)
plot_model(faith1_nodup)
vif(faith1_nodup)
model_faith1_nodup<-plot_model(faith1_nodup)

#use MUMIN::DREDGE to test all models options and find best one and average
faith1_nodup_dredge<-dredge(faith1_nodup)
faith1_nodup_dredge
model.avg(faith1_nodup_dredge)
#plot(rare_faith1_dredge)
averaged_faith1_nodup<-summary(model.avg(get.models(faith1_nodup_dredge, subset = delta < 7)))
averaged_faith1_nodup
sw(averaged_faith1_nodup)


######################### repeat for alleles - using rarefied data########

scatterplotMatrix(~"TLR3","Ase-ua5","Ase-ua3","Ase-ua6","Ase-ua7","Ase-ua9","Ase-ua11","Ase-ua8","Ase-ua4", "Ase-ua1","Ase-ua10", "Ase-dab3","Ase-dab4","Ase-dab5", data=rarealphadata_allC_nodup)

MHC1_CORRELATION<-ggpairs(rarealphadata_allc_nodup,columns=c("TLR3","Ase-ua5","Ase-ua3","Ase-ua6","Ase-ua7","Ase-ua9","Ase-ua11","Ase-ua8","Ase-ua4", "Ase-ua1","Ase-ua10", "Ase-dab3","Ase-dab4","Ase-dab5"))


#Shannon differences - rarefied - alleles ####

rare_shannon1_nodup <- lm(Shannon ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5, data = rarealphadata_allc_nodup, na.action = "na.fail")

#check collinearity for each variable
vif(rare_shannon1_nodup)

#plot residuaals and check effects 
resp = simulateResiduals(rare_shannon1_nodup)
plot(resp, rank = T)
plot_model(rare_shannon1_nodup)

# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = rare_shannon1_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(rare_shannon1_nodup)
Anova(rare_shannon1_nodup)
#visreg(rare_shannon1_nodup)
model_rare_shannon1_nodup<-plot_model(rare_shannon1_nodup)
summ(rare_shannon1_nodup)
vif(rare_shannon1_nodup)

#use MUMIN::DREDGE to test all models options and find best one and average
rare_shannon1_nodup_dredge<-dredge(rare_shannon1_nodup)
rare_shannon1_nodup_dredge
model.avg(rare_shannon1_nodup_dredge)
#plot(shannon1_dredge)
averaged_rare_shannon1_nodup<-summary(model.avg(get.models(rare_shannon1_nodup_dredge, subset = delta < 7)))
averaged_rare_shannon1_nodup
sw(averaged_rare_shannon1_nodup)

#chao1 full model - rarefied - alleles ####

#full model 
rare_chao1a_nodup <- lm(Chao1 ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5, data = rarealphadata_allc_nodup, na.action = "na.fail")

#model not a good fit - try log transforming
rare_chao1_nodup <- lm(logchao ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5, data = rarealphadata_allc_nodup, na.action = "na.fail")

#check collinearity for each variable
vif(rare_chao1_nodup)

#plot residuaals and check effects 
resp = simulateResiduals(rare_chao1_nodup)
plot(resp, rank = T)
plot_model(rare_chao1_nodup)
model_rare_chao1_nodup<-plot_model(rare_chao1_nodup)


# to test if residuals are underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = rare_chao1_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput) 
testOutliers(simulationOutput, alternative = "two.sided") # also tried removing outliers - made no difference to results so kept in
par(mfrow=c(1,1)) 

# general info
summary(rare_chao1_nodup)
Anova(rare_chao1_nodup)
#visreg(rare_chao1_nodup)
plot_model(rare_chao1_nodup)
vif(rare_chao1_nodup)

#use MUMIN::DREDGE to test all models options and find best one and average
rare_chao1_nodup_dredge<-dredge(rare_chao1_nodup)
rare_chao1_nodup_dredge
model.avg(rare_chao1_nodup_dredge)
#plot(chao1_dredge)
averaged_rare_chao1_nodup<-summary(model.avg(get.models(rare_chao1_nodup_dredge, subset = delta < 7)))
averaged_rare_chao1_nodup
sw(averaged_rare_chao1_nodup)

#faith full model - rarefied - alleles ####


#full model 
rare_faith1a_nodup <- lm(PD ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5, data = rarealphadata_allc_nodup, na.action = "na.fail")

#looks like no variance in bird ID - remove
rare_faith1_nodup <- lm(logFaith ~ FieldPeriodID+SexEstimate+AgeClass5+ScaledHs_obs+ TLR3+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5, data = rarealphadata_allc_nodup, na.action = "na.fail")

#check collinearity for each variable
vif(rare_faith1_nodup)

#plot residuaals and check effects - look good
resp = simulateResiduals(rare_faith1_nodup)
plot(resp, rank = T)
plot_model(rare_faith1_nodup)
plot_model(rare_faith1b,type = "re")
resp = simulateResiduals(rare_faith1a_nodup)
plot(resp, rank = T)
model_Rare_faith1<-plot_model(rare_faith1_nodup)


# to test if RESIDUALS ARE underdispersed or overdispersed and examine QQplot- look fine
par(mfrow=c(2,2))
simulationOutput <- simulateResiduals(fittedModel = rare_faith1_nodup, n = 1000)
testDispersion(simulationOutput = simulationOutput, alternative ="less")
testDispersion(simulationOutput = simulationOutput, alternative ="greater")
testUniformity(simulationOutput = simulationOutput)
testOutliers(simulationOutput, alternative = "two.sided")
par(mfrow=c(1,1)) 

# general info
summary(rare_faith1_nodup)
Anova(rare_faith1_nodup)
#visreg(rare_faith1_nodup)
plot_model(rare_faith1_nodup)
vif(rare_faith1_nodup)
model_rare_faith1_nodup<-plot_model(rare_faith1_nodup)

#use MUMIN::DREDGE to test all models options and find best one and average
rare_faith1_nodup_dredge<-dredge(rare_faith1_nodup)
rare_faith1_nodup_dredge
model.avg(rare_faith1_nodup_dredge)
#plot(rare_faith1_dredge)
averaged_rare_faith1_nodup<-summary(model.avg(get.models(rare_faith1_nodup_dredge, subset = delta < 7)))
averaged_rare_faith1_nodup
sw(averaged_rare_faith1_nodup)

#compare rarefied vs none-rarefied models - seem pretty similar so will use rarefied dataset as has fewer factors in model and accounts for sampling depth ####

averaged_rare_shannon1_nodup
averaged_shannon1_nodup

averaged_rare_shannon2_nodup
averaged_shannon2_nodup

averaged_rare_chao2_nodup
averaged_chao2_nodup

averaged_rare_chao1_nodup
averaged_chao1_nodup

averaged_rare_faith1_nodup
averaged_faith1_nodup

averaged_rare_Faith2_nodup
averaged_Faith2_nodup

#compare model outputs
grid.arrange( model_shannon1,model_rare_shannon1,  ncol = 2)
grid.arrange( model_shannon2,model_rare_shannon2,  ncol = 2)
grid.arrange( model_chao1,model_rare_chao1,  ncol = 2)
grid.arrange( model_chao2,model_rare_chao2,  ncol = 2)
grid.arrange( model_faith1,model_Rare_faith1,  ncol = 2)
grid.arrange( model_faith2,model_Rare_faith2,  ncol = 2)


# Fig 3 ####

# Use Estimates and standard errors from linear conditional model-averaged estimates using rarefied data to create Fig 3 and table S6
################################### Beta diversity analysis ########################

# Using cleaned dataset - ie only faecal samples from Cousin - not including chicks, any samples with read depth <10000, no duplicate samples and no samples which are missing MHC data

# Filter taxa ####

# USe Cousin_final_cleaned dataset

# Filter 2: Keep taxa when appearing in minimum % samples
prev0 <- apply(X = otu_table(Cousin_final_cleaned),
               MARGIN = ifelse(taxa_are_rows(Cousin_final_cleaned), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# in 2.5% of samples
Cousin_final_clean_beta <- prune_taxa((prev0 > 5), Cousin_final_cleaned)
ntaxa(Cousin_final_clean_beta) # 1934

# in 2% of samples
Cousin_final_clean_prev2 <- prune_taxa((prev0 > 4), Cousin_final_clean_cleaned)
ntaxa(Cousin_final_clean_prev4) # 2286

# in 1% of samples
Cousin_final_clean_prev1 <- prune_taxa((prev0 > 2), Cousin_final_clean_cleaned)
ntaxa(Cousin_final_clean_prev2) # 3838


#How Many Taxa/sequences did We Lose?    
ntaxa(Cousin_final_cleaned); ntaxa(Cousin_final_clean_beta)
sum(sample_sums(Cousin_final_cleaned)); sum(sample_sums(Cousin_final_clean_beta))
mean(sample_sums(Cousin_final_cleaned)); mean(sample_sums(Cousin_final_clean_beta))


# Figure S2 ####

#Make graph of prevalence of each ASV

# (in how many samples did each ASV appear at least once)

prev0 <- apply(X = otu_table( Cousin_final_mhc_nochick_nodup),
               MARGIN = ifelse(taxa_are_rows(Cousin_final_mhc_nochick_nodup), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts 

prevdf <- data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(Cousin_final_mhc_nochick_nodup),
                     tax_table(Cousin_final_mhc_nochick_nodup))

prevalenceThreshold2 <-2
prevalenceThreshold5 <-4   # 2%
prevalenceThreshold25 <-5    
min_reads10<-10
min_reads100<-100
min_reads50<-50

p_prev <- ggplot(prevdf, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold25, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = min_reads50, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.6) +
  scale_y_log10() + 
  scale_x_log10(labels = c(expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5)), 
                breaks = c(10, 100, 1000, 10000, 100000)) +
  xlab("Total Abundance") +
  facet_wrap(~Phylum) +
    theme_classic() +
  guides(color=FALSE)

p_prev

ggsave("FigS2.jpeg", p_prev, width = 27, height =13 ,units="cm", dpi = 300)


# CSS transformation ####


library("metagenomeSeq")

## Convert the phyloseq object (using 2% filter) to a metagenomeSeq object 
metagenome.obj <- phyloseq_to_metagenomeSeq(Cousin_final_clean_beta)
## Calculate the proper percentile by which to normalize counts
cNstat <- metagenomeSeq::cumNormStatFast(metagenome.obj) 
# cNstat  #0.5
## Normalise counts
metagenome.obj <- metagenomeSeq::cumNorm(metagenome.obj, p = cNstat)
## Export the normalised count table
metag.norm.counts <- metagenomeSeq::MRcounts(metagenome.obj, norm = TRUE)
## Add a pseudocount of 0.0001 to the table and log transform
metag.norm.counts_log <- log(metag.norm.counts+0.0001)
## Substract the value from log of pseudocount to preserve zeros of the original counts
metag.norm.counts_log2 <- metag.norm.counts_log-(log(0.0001))
## Make a new phyloseq object with with the new OTU table
otu_normMG.obj <- otu_table(metag.norm.counts_log2, taxa_are_rows = TRUE)
Cousin_final_clean_beta_css<- phyloseq(otu_normMG.obj, TAX,META,phy_tree)

# repeat for dataset with different prev filters ####

## Convert the phyloseq object to a metagenomeSeq object
metagenome.obj <- phyloseq_to_metagenomeSeq(Cousin_final_clean_prev2)
## Calculate the proper percentile by which to normalize counts
cNstat <- metagenomeSeq::cumNormStatFast(metagenome.obj) 
# cNstat  #0.5
## Normalise counts
metagenome.obj <- metagenomeSeq::cumNorm(metagenome.obj, p = cNstat)
## Export the normalised count table
metag.norm.counts <- metagenomeSeq::MRcounts(metagenome.obj, norm = TRUE)
## Add a pseudocount of 0.0001 to the table and log transform
metag.norm.counts_log <- log(metag.norm.counts+0.0001)
## Substract the value from log of pseudocount to preserve zeros of the original counts
metag.norm.counts_log2 <- metag.norm.counts_log-(log(0.0001))
## Make a new phyloseq object with with the new OTU table
otu_normMG.obj <- otu_table(metag.norm.counts_log2, taxa_are_rows = TRUE)
Cousin_final_clean_prev2_css<- phyloseq(otu_normMG.obj, TAX,META,phy_tree)

## Convert the phyloseq object to a metagenomeSeq object
metagenome.obj <- phyloseq_to_metagenomeSeq(Cousin_final_clean_prev4)
## Calculate the proper percentile by which to normalize counts
cNstat <- metagenomeSeq::cumNormStatFast(metagenome.obj) 
# cNstat  #0.5
## Normalise counts
metagenome.obj <- metagenomeSeq::cumNorm(metagenome.obj, p = cNstat)
## Export the normalised count table
metag.norm.counts <- metagenomeSeq::MRcounts(metagenome.obj, norm = TRUE)
## Add a pseudocount of 0.0001 to the table and log transform
metag.norm.counts_log <- log(metag.norm.counts+0.0001)
## Substract the value from log of pseudocount to preserve zeros of the original counts
metag.norm.counts_log2 <- metag.norm.counts_log-(log(0.0001))
## Make a new phyloseq object with with the new OTU table
otu_normMG.obj <- otu_table(metag.norm.counts_log2, taxa_are_rows = TRUE)
Cousin_final_clean_prev4_css<- phyloseq(otu_normMG.obj, TAX,META,phy_tree)

## Convert the phyloseq object to a metagenomeSeq object
metagenome.obj <- phyloseq_to_metagenomeSeq(Cousin_final_clean_prev5)
## Calculate the proper percentile by which to normalize counts
cNstat <- metagenomeSeq::cumNormStatFast(metagenome.obj) 
# cNstat  #0.5
## Normalise counts
metagenome.obj <- metagenomeSeq::cumNorm(metagenome.obj, p = cNstat)
## Export the normalised count table
metag.norm.counts <- metagenomeSeq::MRcounts(metagenome.obj, norm = TRUE)
## Add a pseudocount of 0.0001 to the table and log transform
metag.norm.counts_log <- log(metag.norm.counts+0.0001)
## Substract the value from log of pseudocount to preserve zeros of the original counts
metag.norm.counts_log2 <- metag.norm.counts_log-(log(0.0001))
## Make a new phyloseq object with with the new OTU table
otu_normMG.obj <- otu_table(metag.norm.counts_log2, taxa_are_rows = TRUE)
Cousin_final_clean_prev5_css<- phyloseq(otu_normMG.obj, TAX,META,phy_tree)

# also rarefy the abundance filtered dataset
set.seed(28367)
Cousin_final_clean_rare = rarefy_even_depth(Cousin_final_clean_beta, sample.size = 10000)

# PERMANOVA =========================================================

# all final permanova models used in manuscript based on variance-stabilised data= Cousin_final_clean_css

#metadata
metadata <- as(sample_data(Cousin_final_clean_beta_css), "data.frame")

# generate distance matrix for each metric
unifrac<- phyloseq::distance(Cousin_final_clean_beta_css, method="unifrac")
weighted_unifrac<- phyloseq::distance(Cousin_final_clean_beta_css, method="wunifrac")

# convert variables to factors
sample_data(Cousin_final_clean_beta_css)$Ase-ua9 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua9)
sample_data(Cousin_final_clean_beta_css)$Ase-ua4 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua4)
sample_data(Cousin_final_clean_beta_css)$Ase-ua5 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua5)
sample_data(Cousin_final_clean_beta_css)$Ase-ua7 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua7)
sample_data(Cousin_final_clean_beta_css)$Ase-ua1 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua1)
sample_data(Cousin_final_clean_beta_css)$Ase-ua3 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua3)
sample_data(Cousin_final_clean_beta_css)$Ase-ua6 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua6)
sample_data(Cousin_final_clean_beta_css)$Ase-ua8 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua8)
sample_data(Cousin_final_clean_beta_css)$Ase-ua11 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-ua11)
sample_data(Cousin_final_clean_beta_css)$Ase-dab3 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-dab3)
sample_data(Cousin_final_clean_beta_css)$Ase-dab4 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-dab4)
sample_data(Cousin_final_clean_beta_css)$Ase-dab5 <- factor(sample_data(Cousin_final_clean_beta_css)$Ase-dab5)

# also tested rarefied data

#metadata
metadata <- as(sample_data(Cousin_final_clean_rare), "data.frame")

# generate distance matrix for each metric
unifrac<- phyloseq::distance(Cousin_final_clean_rare, method="unifrac")
weighted_unifrac<- phyloseq::distance(Cousin_final_clean_rare, method="wunifrac")

# convert variables to factors
sample_data(Cousin_final_clean_rare)$Ase-ua9 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua9)
sample_data(Cousin_final_clean_rare)$Ase-ua4 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua4)
sample_data(Cousin_final_clean_rare)$Ase-ua5 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua5)
sample_data(Cousin_final_clean_rare)$Ase-ua7 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua7)
sample_data(Cousin_final_clean_rare)$Ase-ua1 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua1)
sample_data(Cousin_final_clean_rare)$Ase-ua3 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua3)
sample_data(Cousin_final_clean_rare)$Ase-ua6 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua6)
sample_data(Cousin_final_clean_rare)$Ase-ua8 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua8)
sample_data(Cousin_final_clean_rare)$Ase-ua11 <- factor(sample_data(Cousin_final_clean_rare)$Ase-ua11)
sample_data(Cousin_final_clean_rare)$Ase-dab3 <- factor(sample_data(Cousin_final_clean_rare)$Ase-dab3)
sample_data(Cousin_final_clean_rare)$Ase-dab4 <- factor(sample_data(Cousin_final_clean_rare)$Ase-dab4)
sample_data(Cousin_final_clean_rare)$Ase-dab5 <- factor(sample_data(Cousin_final_clean_rare)$Ase-dab5)


# Run PERMANOVA####

library(vegan)

# set seed so get reproducible results - first try with rarefied data
set.seed(28367)
un_full <- adonis2(unifrac~ AgeClass5+SexEstimate+ TLR3+Hs_obs+FieldPeriodID+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5, data = metadata, by = "margin", perm=9999)
un_full
plot(un_full)

set.seed(28367)
weun_full <- adonis2(weighted_unifrac ~ AgeClass5+SexEstimate+ TLR3+Hs_obs+FieldPeriodID+Ase-ua1+Ase-ua3+Ase-ua4+Ase-ua5+Ase-ua6+Ase-ua7+Ase-ua8+Ase-ua9+Ase-ua11+Ase-dab3+Ase-dab4+Ase-dab5,data = metadata, by = "margin",  perm=9999)
weun_full
plot(weun_full)

# same but for DIVERSITY

set.seed(28367)
un_div <- adonis2(unifrac~ AgeClass5+SexEstimate+ TLR3+Hs_obs+FieldPeriodID+MHC1_Diversity+MHC1_Diversity_2+MHC2_Diversity+MHC2_Diversity_2,data = metadata, by = "margin", perm=9999)
un_div
plot(un_div)

set.seed(28367)
weun_div<- adonis2(weighted_unifrac ~ AgeClass5+SexEstimate+ TLR3+Hs_obs+FieldPeriodID+MHC1_Diversity+MHC1_Diversity_2+MHC2_Diversity+MHC2_Diversity_2, data = metadata, by = "margin", perm=9999)
weun_div
plot(weun_div)

# Are we seeing differences in group means or their dispersion? Use betadisper ####

unifrrac_disp_age <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="unifrac"), group = metadata$AgeClass5)
unifrrac_disp_sex <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="unifrac"), group = metadata$SexEstimate)
unifrrac_disp_Ase-ua4 <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="unifrac"), group = metadata$Ase-ua4)
unifrrac_disp_field <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="unifrac"), group = metadata$FieldPeriodID)
anova(unifrrac_disp_sex)
anova(unifrrac_disp_age)
anova(unifrrac_disp_Ase-ua4)
anova(unifrrac_disp_field)
plot(unifrrac_disp_sex)
plot(unifrrac_disp_field)
boxplot(unifrrac_disp_sex)
boxplot(unifrrac_disp_field)
TukeyHSD(unifrrac_disp_sex)
TukeyHSD(unifrrac_disp_field)
permutest(unifrrac_disp_sex)

weunifrrac_disp_Ase-ua7<- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="weighted_unifrac"), group = metadata$Ase-ua7)
weunifrrac_disp_Ase-ua1 <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="weighted_unifrac"), group = metadata$Ase-ua1)
weunifrrac_disp_Ase-ua11 <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="weighted_unifrac"), group = metadata$Ase-ua11)
weunifrrac_disp_field <- betadisper(d = phyloseq::distance(Cousin_final_clean_css, method="weighed_unifrac"), group = metadata$FieldPeriodID)
anova(weunifrrac_disp_Ase-ua7)
anova(weunifrrac_disp_Ase-ua1)
anova(weunifrrac_disp_Ase-ua11)
anova(weunifrrac_disp_field)
plot(weunifrrac_disp_sex)
plot(weunifrrac_disp_field)
boxplot(weunifrrac_disp_field)
boxplot(weunifrrac_disp_sex)
TukeyHSD(weunifrrac_disp_field)
permutest(weunifrrac_disp_field)

# Figure 4 ####

#plot for unweighted unifran
ord.nmds.UNIFRAC <- ordinate(Cousin_final_clean_beta_css, method="PCoA", distance="unifrac", weighted=FALSE)

#plots for weighted unifrac
ord.nmds.weUNIFRAC <- ordinate(Cousin_final_clean_beta_css, method="PCoA", distance="unifrac", weighted=TRUE)

### ADD COLUMNS TO METADATA

sample_data(Cousin_final_clean_beta_css)$MDS1<-data.frame(ord.nmds.weUNIFRAC$vectors)$Axis.1
sample_data(Cousin_final_clean_beta_css)$MDS2<-data.frame(ord.nmds.weUNIFRAC$vectors)$Axis.2
sample_data(Cousin_final_clean_beta_css)$MDS3<-data.frame(ord.nmds.weUNIFRAC$vectors)$Axis.3
sample_data(Cousin_final_clean_beta_css)$MDS4<-data.frame(ord.nmds.weUNIFRAC$vectors)$Axis.4

# ggplot solution
pca    <-prcomp(df, scale.=T, retx=T)  # principal components analysis
# gg: data frame of PC1 and PC2 scores with corresponding cluster
gg <- data.frame(cluster=factor(km$MHC1.05), x=scores$MDS1, y=scores$MDS2)
# calculate cluster centroid locations
centroids <- aggregate(cbind(x,y)~cluster,data=gg,mean)

# plot cluster map
library(ggplot2)
ggplot(gg, aes(x,y, color=cluster))+
  geom_point(size=3) +
  geom_point(data=centroids, size=4) +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y))+
  geom_path(data=conf.rgn)

ord1<-plot_ordination(Cousin_final_clean_beta_css, ord.nmds.weUNIFRAC, color="Ase-ua7")
PCOA_Aseua7<- ord1 + stat_ellipse(type="norm", size=1) + theme_classic()+ theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16),legend.text = element_text(size=14), legend.title = element_text(size=16))+  theme(legend.position=c(0.8, 0.1))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+scale_color_manual(name="Is Ase-ua7 present?",labels = c("No", "Yes"), values = c("#ffca47","#238A8DFF" ))

ord1<-plot_ordination(Cousin_final_clean_beta_css, ord.nmds.weUNIFRAC, color="Ase-ua11")
PCOA_Aseua11<-ord1 + stat_ellipse(type="norm", size=1) + theme_classic()+ theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16),legend.text = element_text(size=14), legend.title = element_text(size=16))+  theme(legend.position=c(0.8, 0.1))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+scale_color_manual(name="Is Ase_ua11 present?",labels = c("No", "Yes"), values = c("#ffca47","#238A8DFF" ))

ord1<-plot_ordination(Cousin_final_clean_beta_css, ord.nmds.weUNIFRAC, color="Ase-ua1")
PCOA_Aseua1<-ord1 + stat_ellipse(type="norm", size=1) + theme_classic()+ theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16),legend.text = element_text(size=14), legend.title = element_text(size=16))+  theme(legend.position=c(0.8, 0.1))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+scale_color_manual(name="Are Ase_ua1/10 present?",labels = c("No", "Yes"), values = c("#ffca47","#238A8DFF" ))

library(gridExtra)

PCOAplot_Aseua7 <- arrangeGrob(PCOA_Aseua7, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20)))

PCOAplot_Aseua11 <- arrangeGrob(PCOA_Aseua11, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

PCOAplot_Aseua1 <- arrangeGrob(PCOA_Aseua1, top = textGrob("C", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", type ="bold", fontsize=20, fontfamily="Times Roman")))

combo_plot<-grid.arrange( PCOAplot_Aseua7,PCOAplot_Aseua11,PCOAplot_Aseua1, ncol = 3)
ggsave("Figure_4.png", combo_plot, width = 21, height = 7, dpi = 300)

# Figure S5 ####

ord1<-plot_ordination(Cousin_final_clean_beta_css, ord.nmds.weUNIFRAC, color="AgeClass5")
ord1 + stat_ellipse(type="norm", size=1.4) + theme_classic()+ theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16),legend.text = element_text(size=14), legend.title = element_text(size=16))+  theme(legend.position=c(0.88, 0.15))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+scale_color_manual(name="Age class",labels = c("Adult", "Fledgling", "Old fledgling", "Sub-adult"), values = c("#440154FF" , "#FDE725FF","#3CBB75FF" ,"#404788FF"))+
  labs(x="weighted UniFrac Axis 1 (46.9%)", y="weighted UniFrac Axis 2 (10.6%)")

ord1<-plot_ordination(Cousin_final_clean_beta_css, ord.nmds.UNIFRAC, color="AgeClass5")
ord1 + stat_ellipse(type="norm", size=1.4) + theme_classic()+ theme(axis.title.x=element_text(size=20),axis.text.x =element_text(size=16),legend.text = element_text(size=14), legend.title = element_text(size=16))+  theme(legend.position=c(0.15, 1.3))+theme(axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+scale_color_manual(name="Age class",labels = c("Adult", "Fledgling", "Old fledgling", "Sub-adult"), values = c("#440154FF" , "#FDE725FF","#3CBB75FF" ,"#404788FF"))+
  labs(x="unweighted UniFrac Axis 1 (11.9%)", y="unweighted UniFrac Axis 2 (6.1%)")


########################### De-seq Analysis  #################################
# Variance stabilizing transformation ==========================================

library(DESeq2)

# convert to deseq - using filtered dataset from beta diversity analysis
nes_dds <- phyloseq_to_deseq2(Cousin_final_clean_beta, ~1)
# estimate size factors with not including 0 in geometric mean calc
nes_dds  <- estimateSizeFactors(nes_dds, type = "poscounts") %>% 
  estimateDispersions(fitType = "local")
# create new phyloseq object with variance stabilised ASV table
Cousin_final_clean_vst <- Cousin_final_clean_beta
otu_table(Cousin_final_clean_vst) <- otu_table(getVarianceStabilizedData(nes_dds), 
                                               taxa_are_rows = TRUE)


# first remove taxa which have no count - but should already be removed
any(taxa_sums(Cousin_final_clean_beta) == 0)
Cousin_final_clean_beta <- prune_taxa(taxa_sums(Cousin_final_clean_beta) > 0, Cousin_final_clean_beta)
any(taxa_sums(Cousin_final_clean_beta) == 0)

# change presence/absence to a factor rather than continuos variable as coded as 0/1
sample_data(Cousin_final_clean_beta)$Ase-ua7 <- factor(sample_data(Cousin_final_clean_beta)$Ase-ua7)
sample_data(Cousin_final_clean_beta)$Ase-ua11 <- factor(sample_data(Cousin_final_clean_beta)$Ase-ua11)
sample_data(Cousin_final_clean_beta)$Ase-ua1 <- factor(sample_data(Cousin_final_clean_beta)$Ase-ua1)


# De-seq ####

# create deseq object using cleaned data - only included factors which were significant in beta diversity analysis

dsAll <- phyloseq_to_deseq2(Cousin_final_clean_beta, ~ FieldPeriodID+AgeClass5+Ase-ua7+Ase-ua11+Ase-ua1)
# estimate size factors with not including 0 in geometric mean calc
dsAll_vst  <- estimateSizeFactors(dsAll, type = "poscounts")
# deseq analysis - have used local rather than parametric
dsAll_test<- DESeq(dsAll_vst, test="Wald", fitType="local")

# not all genes converged - see message below
#5 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#-- replacing outliers and refitting for 1003 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
#1 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

# see which genes did not converge and find ASV name
View(mcols(dsAll_test)$betaConv)
View(mcols(dsAll_test)@rownames)

# to remove not converged taxa
ddsClean <- dsAll_test[which(mcols(dsAll_test)$betaConv),]

#to plot specific ASVs e.g., 
plotCounts(ddsClean, gene="3660d74cf74458cc9dd7a5e0e075752f", intgroup="Ase-ua11", normalized = F)

# can try to increase iterations used to get full convergence
#dds <- estimateDispersions(dsAll_vst, fitType="local")
#dsAll_test <- nbinomWaldTest(dds, maxit=10000)

#check how well model worked etc
plotDispEsts(dsAll_test) 
plotDispEsts(ddsClean) 

# MHC data####

# Ase-a11
res_Ase-ua11<- results(ddsClean,contrast=c("Ase-ua11","0","1"))
#how many with p<0.01
sum( res_Ase-ua11$padj < 0.01, na.rm=TRUE )
plotMA( res_Ase-ua11 )

# results table
alpha <- 0.01
sigtab_Ase-ua11 <- res_Ase-ua11[which(res_Ase-ua11$padj < alpha), ]
sigtab_Ase-ua11<- cbind(as(sigtab_Ase-ua11, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Ase-ua11), ], "matrix"))
head(sigtab_Ase-ua11)
dim(sigtab_Ase-ua11)
# 44 taxa

## Order results by the largest fold change
x_MHC <- tapply(sigtab_Ase-ua11$log2FoldChange, sigtab_Ase-ua11$Phylum, function(x_MHC) max(x_MHC))
x_MHC <- sort(x_MHC, TRUE)
sigtab_Ase-ua11$Phylum <- factor(as.character(sigtab_Ase-ua11$Phylum), levels=names(x_MHC))

paste( "Overall, we find", length(sigtab_Ase-ua11$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Ase-ua11$log2FoldChange < 0)), "being significantly more abundant when Ase-ua11 allele is absent and", length(which(sigtab_Ase-ua11$log2FoldChange > 0)), "when present.")

write.table(sigtab_Ase-ua11, "DeSeq_results_Ase-UA11.txt")

# Ase-ua7
res_Ase-ua7<- results(ddsClean,contrast=c("Ase-ua7","0","1"))
#how many with p<0.01
sum( res_Ase-ua7$padj < 0.01, na.rm=TRUE )
plotMA( res_Ase-ua7 )

# results table
alpha <- 0.01
sigtab_Ase-ua7 <- res_Ase-ua7[which(res_Ase-ua7$padj < alpha), ]
sigtab_Ase-ua7<- cbind(as(sigtab_Ase-ua7, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Ase-ua7), ], "matrix"))
head(sigtab_Ase-ua7)
dim(sigtab_Ase-ua7)
# 52 taxa


## Order results by the largest fold change
x_MHC <- tapply(sigtab_Ase-ua7$log2FoldChange, sigtab_Ase-ua7$Phylum, function(x_MHC) max(x_MHC))
x_MHC <- sort(x_MHC, TRUE)
sigtab_Ase-ua75$Phylum <- factor(as.character(sigtab_Ase-ua7$Phylum), levels=names(x_MHC))

paste( "Overall, we find", length(sigtab_Ase-ua7$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Ase-ua7$log2FoldChange < 0)), "being significantly more abundant when Ase-ua7 allele is absent and", length(which(sigtab_Ase-ua7$log2FoldChange > 0)), "when present.")

write.table(sigtab_MHC1.05, "DeSeq_results_Ase_ua7.txt")

#Ase-ua1/10
res_Ase-ua1<- results(ddsClean,contrast=c("Ase-ua1","0","1"))
#how many with p<0.01
sum( res_Ase-ua1$padj < 0.01, na.rm=TRUE )
plotMA( res_Ase-ua1 )

# results table
alpha <- 0.01
sigtab_Ase-ua1 <- res_Ase-ua1[which(res_Ase-ua1$padj < alpha), ]
sigtab_Ase-ua1<- cbind(as(sigtab_Ase-ua1, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Ase-ua1), ], "matrix"))
head(sigtab_Ase-ua1)
dim(sigtab_Ase-ua1)
# 99 taxa


## Order results by the largest fold change
x_MHC <- tapply(sigtab_Ase-ua1$log2FoldChange, sigtab_Ase-ua1$Phylum, function(x_MHC) max(x_MHC))
x_MHC <- sort(x_MHC, TRUE)
sigtab_Ase-ua1$Phylum <- factor(as.character(sigtab_Ase-ua1$Phylum), levels=names(x_MHC))

paste( "Overall, we find", length(sigtab_Ase-ua1$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Ase-ua1$log2FoldChange < 0)), "being significantly more abundant when Ase-ua1/10 allele is absent and", length(which(sigtab_Ase-ua1$log2FoldChange > 0)), "when present.")

write.table(sigtab_Ase-ua1, "DeSeq_results_Ase_1.txt")

# Fig 5 ####

#to get class names
View(merged.top.phylum@tax_table@.Data)

## Assign colours to the phyla (matching those from the relative abundance plot)
phylcols <- c(Proteobacteria = "dodgerblue",Firmicutes = "firebrick1",Actinobacteria = "gold",Chloroflexi = "darkorchid",Bacteroidetes = "chartreuse4",Planctomycetes = "aquamarine",Patescibacteria = "darkorange1", Acidobacteria = "darkslategrey",  Verrucomicrobia = "deeppink1", Cyanobacteria = "greenyellow", Deinococcus_Thermus = "grey80",Notassigned="grey80",  Fusobacteria = "grey80", Rokubacteria = "grey80")
## Change name of Candidatus_Saccharibacteria and Deinococcus_Thermus back to the original names used in the table

names(phylcols)[which(names(phylcols)=="Deinococcus_Thermus")] <- "Deinococcus-Thermus"

#export tables then
#import back into R - with missing orders added in for y-axis consistency between plots
plot_Ase_ua11 <- read_csv("DeSeq_results_Ase_ua11_reordered.csv")
plot_Ase_ua7 <- read_csv("DeSeq_results_Ase_ua7_reordered.csv")
plot_Ase_ua1 <- read_csv("DeSeq_results_Ase_1_reordered.csv")

## Make the plot

plot1<-plot_Ase_ua11 %>%  mutate(order=ordered(order,levels=(c("Not assigned","Fusobacteriales","Deinococcales","Chthoniobacterales","Saccharimonadales","Bacteroidales","Uncultured planctomycete","Planctomycetales","Pirellulales","Isosphaerales","Gemmatales","Thermomicrobiales","Chloroflexales","Streptomycetales","Solirubrobacterales","Rubrobacterales","Pseudonocardiales","Propionibacteriales","Microtrichales","Micromonosporales","Micrococcales","Kineosporiales","Frankiales","Corynebacteriales","Lactobacillales","Erysipelotrichales","Clostridiales","Bacillales","Xanthomonadales","Uncultured proteobacteria","Sphingomonadales","Rhodobacterales","Rhizobiales","Reyranellales","Pseudomonadales","Oceanospirillales","Enterobacteriales","Caulobacterales","Acetobacterales")))) %>% arrange(order) %>%
  ggplot( aes(x=log2FoldChange, y=order, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("Ase-ua11 ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot2<-plot_Ase_ua7 %>%  mutate(order=ordered(order,levels=(c("Not assigned","Fusobacteriales","Deinococcales","Chthoniobacterales","Saccharimonadales","Bacteroidales","Uncultured planctomycete","Planctomycetales","Pirellulales","Isosphaerales","Gemmatales","Thermomicrobiales","Chloroflexales","Streptomycetales","Solirubrobacterales","Rubrobacterales","Pseudonocardiales","Propionibacteriales","Microtrichales","Micromonosporales","Micrococcales","Kineosporiales","Frankiales","Corynebacteriales","Lactobacillales","Erysipelotrichales","Clostridiales","Bacillales","Xanthomonadales","Uncultured proteobacteria","Sphingomonadales","Rhodobacterales","Rhizobiales","Reyranellales","Pseudomonadales","Oceanospirillales","Enterobacteriales","Caulobacterales","Acetobacterales")))) %>% arrange(order) %>%
  ggplot( aes(x=log2FoldChange, y=order, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14),  axis.text.y = element_blank())+theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("Ase-ua7")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot3<-plot_Ase_ua1 %>%  mutate(order=ordered(order,levels=(c("Not assigned","","Fusobacteriales","Deinococcales","Chthoniobacterales","Saccharimonadales","Bacteroidales","Uncultured planctomycete","Planctomycetales","Pirellulales","Isosphaerales","Gemmatales","Thermomicrobiales","Chloroflexales","Streptomycetales","Solirubrobacterales","Rubrobacterales","Pseudonocardiales","Propionibacteriales","Microtrichales","Micromonosporales","Micrococcales","Kineosporiales","Frankiales","Corynebacteriales","Lactobacillales","Erysipelotrichales","Clostridiales","Bacillales","Xanthomonadales","Uncultured proteobacteria","Sphingomonadales","Rhodobacterales","Rhizobiales","Reyranellales","Pseudomonadales","Oceanospirillales","Enterobacteriales","Caulobacterales","Acetobacterales")))) %>% arrange(order) %>%
  ggplot( aes(x=log2FoldChange, y=order, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+ 
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_blank())+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  ggtitle("Ase-ua1/10 ")+theme(legend.position = "none")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# to get legend
correct.order <- c("Proteobacteria","Firmicutes","Actinobacteria", "Chloroflexi","Planctomycetes","Bacteroidetes","Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Deinococcus-Thermus","Fusobacteria")
plot_Ase_ua1$Phylum <- factor(plot_Ase_ua1$Phylum,levels = correct.order)

legend_plot<-plot_Ase_ua1 %>%  mutate(Pylum=ordered(Phylum,levels=(c("Proteobacteria","Firmicutes","Actinobacteria", "Chloroflexi","Planctomycetes","Bacteroidetes","Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Deinococcus-Thermus","Fusobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=order, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  ggtitle("Ase-ua1/10 ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(cowplot)
legend <- cowplot::get_legend(legend_plot)

grid.newpage()
grid.draw(legend)

combo_plot<-grid.arrange( plot1,plot2,plot3,legend_plot, ncol = 4)
ggsave("Figure_4.png", combo_plot, width = 21, height = 7, dpi = 300)

# To extract mean relative abundance ####

relabund <- transform_sample_counts(Cousin_final_clean_beta, function(x){(x / sum(x))*100})
dat <- psmelt(relabund)
Ase_ua7_aggregate<-aggregate(Abundance ~ OTU+ Ase-ua7, data=dat, mean)
Ase_ua11_aggregate<-aggregate(Abundance ~ OTU+ Ase-ua11, data=dat, mean)
Ase_ua1_aggregate<-aggregate(Abundance ~ OTU+ Ase-ua1, data=dat, mean)

write.table(Ase_ua7_aggregate, "Ase_ua7_aggregate.txt")
write.table(Ase_ua11_aggregate, "Ase_ua11_aggregate.txt")
write.table(Ase_ua1_aggregate, "Ase_ua1_aggregate.txt")

# Add relative abundance values to appropriate tabs in additional file 2

# Age class  data####

# fl VS ofl
res_Age_OFL_FL<- results(ddsClean,contrast=c("AgeClass5","FL","OFL"))
#how many with p<0.01
sum( res_Age_OFL_FL$padj < 0.01, na.rm=TRUE )
plotMA( res_Age_OFL_FL )

# results table
alpha <- 0.01
sigtab_Age_FL_OFL <- res_Age_OFL_FL[which(res_Age_OFL_FL$padj < alpha), ]
sigtab_Age_FL_OFL<- cbind(as(sigtab_Age_FL_OFL, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Age_FL_OFL), ], "matrix"))
head(sigtab_Age_FL_OFL)
dim(sigtab_Age_FL_OFL)
# 174 taxa
#View(sigtab_Age_FL_OFL)

## Order results by the largest fold change
x_age <- tapply(sigtab_Age_FL_OFL$log2FoldChange, sigtab_Age_FL_OFL$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_Age_FL_OFL$Phylum <- factor(as.character(sigtab_Age_FL_OFL$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_Age_FL_OFL$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Age_FL_OFL$log2FoldChange < 0)), "being significantly more abundant in fledlging and", length(which(sigtab_Age_FL_OFL$log2FoldChange > 0)), "in old fledglings.")

write.table(sigtab_Age_FL_OFL, "DeSeq_results_age_ofl_fl.txt")

#SA VS ofl
res_Age_SA_OFL<- results(ddsClean,contrast=c("AgeClass5","SA","OFL"))
#how many with p<0.01
sum( res_Age_SA_OFL$padj < 0.01, na.rm=TRUE )
plotMA( res_Age_SA_OFL )

# results table
alpha <- 0.01
sigtab_Age_SA_OFL <- res_Age_SA_OFL[which(res_Age_SA_OFL$padj < alpha), ]
sigtab_Age_SA_OFL<- cbind(as(sigtab_Age_SA_OFL, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Age_SA_OFL), ], "matrix"))
head(sigtab_Age_SA_OFL)
dim(sigtab_Age_SA_OFL)
# 193 taxa
#View(sigtab_Age_SA_OFL)

## Order results by the largest fold change
x_age <- tapply(sigtab_Age_SA_OFL$log2FoldChange, sigtab_Age_SA_OFL$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_Age_SA_OFL$Phylum <- factor(as.character(sigtab_Age_SA_OFL$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_Age_SA_OFL$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Age_SA_OFL$log2FoldChange < 0)), "being significantly more abundant in SA and", length(which(sigtab_Age_SA_OFL$log2FoldChange > 0)), "in old fledglings.")

write.table(sigtab_Age_SA_OFL, "DeSeq_results_age_ofl_SA.txt")

#A VS ofl
res_Age_A_OFL<- results(ddsClean,contrast=c("AgeClass5","A","OFL"))
#how many with p<0.01
sum( res_Age_A_OFL$padj < 0.01, na.rm=TRUE )
plotMA( res_Age_A_OFL )

# results table
alpha <- 0.01
sigtab_Age_A_OFL <- res_Age_A_OFL[which(res_Age_A_OFL$padj < alpha), ]
sigtab_Age_A_OFL<- cbind(as(sigtab_Age_A_OFL, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Age_A_OFL), ], "matrix"))
head(sigtab_Age_A_OFL)
dim(sigtab_Age_A_OFL)
# 186 taxa
#View(sigtab_Age_A_OFL)

## Order results by the largest fold change
x_age <- tapply(sigtab_Age_A_OFL$log2FoldChange, sigtab_Age_A_OFL$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_Age_A_OFL$Phylum <- factor(as.character(sigtab_Age_A_OFL$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_Age_A_OFL$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Age_A_OFL$log2FoldChange < 0)), "being significantly more abundant in A and", length(which(sigtab_Age_A_OFL$log2FoldChange > 0)), "in old fledglings.")

write.table(sigtab_Age_A_OFL, "DeSeq_results_age_ofl_A.txt")

#A VS SA
res_Age_A_SA<- results(ddsClean,contrast=c("AgeClass5","A","SA"))
#how many with p<0.01
sum( res_Age_A_SA$padj < 0.01, na.rm=TRUE )
plotMA( res_Age_A_SA )

# results table
alpha <- 0.01
sigtab_Age_A_SA <- res_Age_A_SA[which(res_Age_A_SA$padj < alpha), ]
sigtab_Age_A_SA<- cbind(as(sigtab_Age_A_SA, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Age_A_SA), ], "matrix"))
head(sigtab_Age_A_SA)
dim(sigtab_Age_A_SA)
# 108 taxa
#View(sigtab_Age_A_SA)

## Order results by the largest fold change
x_age <- tapply(sigtab_Age_A_SA$log2FoldChange, sigtab_Age_A_SA$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_Age_A_SA$Phylum <- factor(as.character(sigtab_Age_A_SA$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_Age_A_SA$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Age_A_SA$log2FoldChange < 0)), "being significantly more abundant in A and", length(which(sigtab_Age_A_SA$log2FoldChange > 0)), "in SA.")

write.table(sigtab_Age_A_SA, "DeSeq_results_age_SA_A.txt")

#A VS FL
res_Age_A_FL<- results(ddsClean,contrast=c("AgeClass5","A","FL"))
#how many with p<0.01
sum( res_Age_A_FL$padj < 0.01, na.rm=TRUE )
plotMA( res_Age_A_FL )

# results table
alpha <- 0.01
sigtab_Age_A_FL <- res_Age_A_FL[which(res_Age_A_FL$padj < alpha), ]
sigtab_Age_A_FL<- cbind(as(sigtab_Age_A_FL, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Age_A_FL), ], "matrix"))
head(sigtab_Age_A_FL)
dim(sigtab_Age_A_FL)
# 55 taxa
#View(sigtab_Age_A_FL)

## Order results by the largest fold change
x_age <- tapply(sigtab_Age_A_FL$log2FoldChange, sigtab_Age_A_FL$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_Age_A_FL$Phylum <- factor(as.character(sigtab_Age_A_FL$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_Age_A_FL$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Age_A_FL$log2FoldChange < 0)), "being significantly more abundant in A and", length(which(sigtab_Age_A_FL$log2FoldChange > 0)), "in FL.")

write.table(sigtab_Age_A_FL, "DeSeq_results_age_FL_A.txt")

#SA VS FL
res_Age_SA_FL<- results(ddsClean,contrast=c("AgeClass5","SA","FL"))
#how many with p<0.01
sum( res_Age_SA_FL$padj < 0.01, na.rm=TRUE )
plotMA( res_Age_SA_FL )

# results table
alpha <- 0.01
sigtab_Age_SA_FL <- res_Age_A_FL[which(res_Age_SA_FL$padj < alpha), ]
sigtab_Age_SA_FL<- cbind(as(sigtab_Age_SA_FL, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_Age_SA_FL), ], "matrix"))
head(sigtab_Age_SA_FL)
dim(sigtab_Age_SA_FL)
# 102 taxa
#View(sigtab_Age_SA_FL)

## Order results by the largest fold change
x_age <- tapply(sigtab_Age_SA_FL$log2FoldChange, sigtab_Age_SA_FL$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_Age_SA_FL$Phylum <- factor(as.character(sigtab_Age_SA_FL$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_Age_SA_FL$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_Age_SA_FL$log2FoldChange < 0)), "being significantly more abundant in SA and", length(which(sigtab_Age_SA_FL$log2FoldChange > 0)), "in FL.")

write.table(sigtab_Age_SA_FL, "DeSeq_results_age_FL_SA.txt")

# Figure S6 ####

#import back into R - field season
plot_OFL_FL <- read_csv("DeSeq_results_age_ofl_fl_reordered.csv")
plot_OFL_SA <- read_csv("DeSeq_results_age_ofl_SA_reordered.csv")
plot_OFL_A <- read_csv("DeSeq_results_age_ofl_A_reordered.csv")
plot_FL_SA <- read_csv("DeSeq_results_age_FL_SA_reordered.csv")
plot_FL_A <- read_csv("DeSeq_results_age_FL_A_REORDERED.csv")
plot_SA_A <- read_csv("DeSeq_results_age_SA_A_REORDERED.csv")

## Make the plot

plot1<-plot_OFL_FL %>% mutate(Phylum=ordered(Phylum,levels=(c("Rokubacteria","Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("FL -  OFL ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plot2<-plot_OFL_SA %>% mutate(Phylum=ordered(Phylum,levels=(c("Rokubacteria","Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("SA - OFL")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot3<-plot_OFL_A %>% mutate(Phylum=ordered(Phylum,levels=(c("Rokubacteria","Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("A - OFL")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
## Make the plot

plot4<-plot_FL_SA %>% mutate(Phylum=ordered(Phylum,levels=(c("Rokubacteria","Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("FL -  SA ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot5<-plot_FL_A %>% mutate(Phylum=ordered(Phylum,levels=(c("Rokubacteria","Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("FL - A")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot6<-plot_SA_A %>% mutate(Phylum=ordered(Phylum,levels=(c("Rokubacteria","Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("SA - A")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# to get legend
correct.order <- c("Proteobacteria","Firmicutes","Actinobacteria", "Chloroflexi","Planctomycetes","Bacteroidetes","Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Deinococcus-Thermus","Fusobacteria")
plot_plot_OFL_FL$Phylum <- factor(plot_plot_OFL_FL$Phylum,levels = correct.order)

legend_plot<-plot_plot_OFL_FL %>%  mutate(Pylum=ordered(Phylum,levels=(c("Proteobacteria","Firmicutes","Actinobacteria", "Chloroflexi","Planctomycetes","Bacteroidetes","Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Deinococcus-Thermus","Fusobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=order, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  ggtitle("Ase-ua1/10 ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(cowplot)
legend <- cowplot::get_legend(legend_plot)

grid.newpage()
grid.draw(legend)

combo_plot<-grid.arrange( legend_plot,plot1,plot2,legend_plot,plot3,plot4,legend_plot,plot5,plot6, ncol = 3)
ggsave("Figure_S6.png", combo_plot, width = 21, height = 7, dpi = 300)

# Field period ####

#Summer 2017 and 18
res_FP_S2017_S2018<- results(ddsClean,contrast=c("FieldPeriodID","Summer2017","Summer2018"))
#how many with p<0.01
sum(res_FP_S2017_S2018$padj < 0.01, na.rm=TRUE )
plotMA( res_FP_S2017_S2018 )

# results table
alpha <- 0.01
sigtab_FP_S2017_S2018 <- res_FP_S2017_S2018[which(res_FP_S2017_S2018$padj < alpha), ]
sigtab_FP_S2017_S2018<- cbind(as(sigtab_FP_S2017_S2018, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_FP_S2017_S2018), ], "matrix"))
head(sigtab_FP_S2017_S2018)
dim(sigtab_FP_S2017_S2018)
# 144 taxa
#View(sigtab_FP_S2017_S2018)

## Order results by the largest fold change
x_field <- tapply(sigtab_FP_S2017_S2018$log2FoldChange, sigtab_FP_S2017_S2018$Phylum, function(x_field) max(x_field))
x_field <- sort(x_field, TRUE)
sigtab_FP_S2017_S2018$Phylum <- factor(as.character(sigtab_FP_S2017_S2018$Phylum), levels=names(x_field))

paste( "Overall, we find", length(sigtab_FP_S2017_S2018$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_FP_S2017_S2018$log2FoldChange < 0)), "being significantly more abundant in Summer2017 and", length(which(sigtab_FP_S2017_S2018$log2FoldChange > 0)), "in Summer 2018.")

write.table(sigtab_FP_S2017_S2018, "DeSeq_results_FP_S2017_S2018.txt")

#Summer 2018 and winter 18
res_FP_S2018_w2018<- results(ddsClean,contrast=c("FieldPeriodID","Summer2018","Winter2018"))
#how many with p<0.01
sum(res_FP_S2018_w2018$padj < 0.01, na.rm=TRUE )
plotMA( res_FP_S2018_w2018 )

# results table
alpha <- 0.01
sigtab_FP_S2018_w2018 <- res_FP_S2018_w2018[which(res_FP_S2018_w2018$padj < alpha), ]
sigtab_FP_S2018_w2018<- cbind(as(sigtab_FP_S2018_w2018, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_FP_S2018_w2018), ], "matrix"))
head(sigtab_FP_S2018_w2018)
dim(sigtab_FP_S2018_w2018)
# 225 taxa
#View(sigtab_FP_S2018_w2018)

## Order results by the largest fold change
x_age <- tapply(sigtab_FP_S2018_w2018$log2FoldChange, sigtab_FP_S2018_w2018$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_FP_S2018_w2018$Phylum <- factor(as.character(sigtab_FP_S2018_w2018$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_FP_S2018_w2018$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_FP_S2018_w2018$log2FoldChange < 0)), "being significantly more abundant in Summer2018 and", length(which(sigtab_FP_S2018_w2018$log2FoldChange > 0)), "in Winter 2018.")

write.table(sigtab_FP_S2018_w2018, "DeSeq_results_FP_S2018_w2018.txt")

#Summer 2017 and winter 18
res_FP_S2017_w2018<- results(ddsClean,contrast=c("FieldPeriodID","Summer2017","Winter2018"))
#how many with p<0.01
sum(res_FP_S2017_w2018$padj < 0.01, na.rm=TRUE )
plotMA( res_FP_S2017_w2018 )

# results table
alpha <- 0.01
sigtab_FP_S2017_w2018 <- res_FP_S2017_w2018[which(res_FP_S2017_w2018$padj < alpha), ]
sigtab_FP_S2017_w2018<- cbind(as(sigtab_FP_S2017_w2018, "data.frame"), as(tax_table(Cousin_final_clean_beta)[rownames(sigtab_FP_S2017_w2018), ], "matrix"))
head(sigtab_FP_S2017_w2018)
dim(sigtab_FP_S2017_w2018)
# 227 taxa
#View(sigtab_FP_S2017_w2018)

## Order results by the largest fold change
x_age <- tapply(sigtab_FP_S2017_w2018$log2FoldChange, sigtab_FP_S2017_w2018$Phylum, function(x_age) max(x_age))
x_age <- sort(x_age, TRUE)
sigtab_FP_S2017_w2018$Phylum <- factor(as.character(sigtab_FP_S2017_w2018$Phylum), levels=names(x_age))

paste( "Overall, we find", length(sigtab_FP_S2017_w2018$log2FoldChange)," significantly differentially abundant OTUs with", length(which(sigtab_FP_S2017_w2018$log2FoldChange < 0)), "being significantly more abundant in Summer2017 and", length(which(sigtab_FP_S2017_w2018$log2FoldChange > 0)), "in Winter 2018.")

write.table(sigtab_FP_S2017_w2018, "DeSeq_results_FP_S2017_w2018.txt")

# Fig S7####

#import back into R - field season
plot_S2017_W2018 <- read_csv("DeSeq_results_FP_S2017_w2018_reordered.csv")
plot_S2018_W2018 <- read_csv("DeSeq_results_FP_S2018_w2018_reordered.csv")
plot_S2017_S2018 <- read_csv("DeSeq_results_FP_S2017_S2018_reordered.csv")

## Make the plot

plot1<-plot_S2017_W2018 %>% mutate(Phylum=ordered(Phylum,levels=(c("Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("Major 2017 -  Minor 2018 ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot2<-plot_S2017_S2018 %>% mutate(Phylum=ordered(Phylum,levels=(c("Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("Major 2017 - Major 2018")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot3<-plot_S2018_W2018 %>% mutate(Phylum=ordered(Phylum,levels=(c("Fusobacteria","Deinococcus-Thermus","Cyanobacteria","Verrucomicrobia","Acidobacteria","Patescibacteria","Bacteroidetes","Planctomycetes","Chloroflexi","Actinobacteria","Firmicutes","Proteobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=Phylum, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  theme(legend.text = element_text( size = 10),legend.title = element_text(face="bold"))+
  ggtitle("Major 2018 - Minor 2018")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# to get legend
correct.order <- c("Proteobacteria","Firmicutes","Actinobacteria", "Chloroflexi","Planctomycetes","Bacteroidetes","Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Deinococcus-Thermus","Fusobacteria")
plot_S2017_W2018$Phylum <- factor(plot_S2017_W2018$Phylum,levels = correct.order)

legend_plot<-plot_S2017_W2018 %>%  mutate(Pylum=ordered(Phylum,levels=(c("Proteobacteria","Firmicutes","Actinobacteria", "Chloroflexi","Planctomycetes","Bacteroidetes","Patescibacteria","Acidobacteria","Verrucomicrobia","Cyanobacteria","Deinococcus-Thermus","Fusobacteria")))) %>% arrange(Phylum) %>%
  ggplot( aes(x=log2FoldChange, y=order, colour=Phylum)) +
  geom_point(size=2.5) + 
  geom_vline(xintercept = 0,linetype = 2, colour="gray44")+
  theme_bw()+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 5.5, linetype=0), ncol = 1))+
  ggtitle("Ase-ua1/10 ")+
  scale_colour_manual(values=phylcols)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(cowplot)
legend <- cowplot::get_legend(legend_plot)

grid.newpage()
grid.draw(legend)

combo_plot<-grid.arrange( plot1,plot2,plot3,legend_plot, ncol = 4)
ggsave("Figure_S7.png", combo_plot, width = 21, height = 7, dpi = 300)
