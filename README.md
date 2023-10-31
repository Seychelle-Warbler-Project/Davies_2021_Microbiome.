# Davies_2021_Microbiome

Data and script needed to reproduce the results in the paper Davies et al (2022): Immunogenetic variation shapes the gut microbiome in a natural vertebrate population, doi:10.1186/s40168-022-01233-y

All variable names in these data should be self-explanatory, but please don't hesitate to send an email to the corresponding author if anything is unclear

# R script files:
## R_script_Davies_et_al_GM
This R script contains all the R code needed to reproduce the results from the paper - using the data from the four files in the folder - Phyloseq_files

## DECONTAM_R_script_Davies_et_al_GM
This R script contains  R code needed to reproduce the decontamination steps from the paper - using the data from the files found in the folder: Contamination_included


# Folders included:

## Phyloseq_files
This folder contains the decontaminated raw data files needed for R_script_Davies_et_al_GM. Microbiome files are exported directly from QIIME2 and MHC data generated using AMPLISAS

merged_otu_table.txt: This file contains the raw abundances of Amplicon Sequence Variants (ASV) inferred for each faecal sample using QIIMEII analysis. The paired-end sequences for each sample used to generate the table have been deposited with the ENA under the accession numbers PRJEB45408. The ASV table is used to generate phyloseq objects

taxonomy_tsv: This file contains the taxonomy of ASVs in ASV_table.csv. It was generated as part of the QIIMEII pipeline. This file is used to generate phyloseq objects

tree.nwk: A rooted phylogenetic tree of the ASVs in ASV_table.csv. It was generated as part of the QIIMEII pipeline. Used to generate phyloseq objects and conduct phylogenetically aware analyses

Sample_metadata_Daviesetal.csv: This file contains the sample metadata - All variable names in these data should be self-explanatory, but please don't hesitate to send an email to the corresponding author if anything is unclear.

## Phyloseq_files_Contamination_included
This folder contains the data files needed for the decontamination script 

## Tidied_tables_for_plots
This folder contains files which have been tidied for aesthetic purposes and are used in the R_script_Davies_et_al_GM  


