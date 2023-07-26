# Check for Kat & Sean 
library(tidyverse)

# Song et al. 
ASV_S = read.table("Data/ForSean/Song_ASVs_counts_filt.txt")
Taxo_S = read.table("Data/ForSean/Song_ASVs_taxonomy_filt.txt")
Meta_S = read.table("Data/ForSean/Song_Metadata.txt")

dim(ASV_S)
dim(Taxo_S)
dim(Meta_S)

setdiff(colnames(ASV_S),Meta_S$SampleID) # Sample IDs correspond in both tables 
setdiff(rownames(ASV_S),Taxo_S$ASV) # ASVs corresponds in both tables 


# Youngblut et al. 
ASV_Y = read.table( "Data/ForSean/Youngblut_ASVs_counts_filt.txt")
Taxo_Y = read.table( "Data/ForSean/Youngblut_ASVs_taxonomy_filt.txt")
Meta_Y = read.table( "Data/ForSean/Youngblut_Metadata.txt")


dim(ASV_Y)
dim(Taxo_Y)
dim(Meta_Y)

setdiff(colnames(ASV_Y),Meta_Y$Sample_ID) # Sample IDs correspond in both tables 
setdiff(rownames(ASV_Y),Taxo_Y$ASV) # ASVs corresponds in both tables 
