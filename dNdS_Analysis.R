# Title     : TODO
# Objective : TODO
# Created by: schencro
# Created on: 12/18/17

library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")
library(ggplot2)

data("dataset_normalskin_genes", package="dndscv")
data("dataset_normalskin", package="dndscv")
m <- subset(m, m$alt %in% c("A","C","G","T"))
m <- subset(m, m$ref %in% c("A","C","G","T"))

#dndsskin_2r = dndscv(m, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, sm = "2r_3w")
dndsskin = dndscv(m, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
sel_cv = dndsskin$sel_cv
print(head(sel_cv[sel_cv$qglobal_cv<0.1, c(1:10,19)]), digits = 3)
print(dndsskin$globaldnds, digits = 5)

#------Patient Data------#
load("~/Desktop/Darryl_collab/Patient_Data/Patient_Data_Shiny/For_Shiny_App/Patient_Data/df_patients.RData")
imptPat <- data.frame(
sampleID=as.character(df_patients$sample),
chr=as.character(df_patients$chr),
pos=as.integer(df_patients$pos),
ref=as.character(df_patients$ref_nt),
mut=as.character(df_patients$mut_nt)
)
colnames(imptPat) <- c("sampleID","chr","pos","ref","alt"  )
imptPat <- subset(imptPat, imptPat$alt %in% c("A","C","G","T"))
imptPat <- subset(imptPat, imptPat$ref %in% c("A","C","G","T"))
imptPat$sampleID <- as.character(imptPat$sampleID)
imptPat$chr <- as.character(imptPat$chr)
imptPat$ref <- as.character(imptPat$ref)
imptPat$alt <- as.character(imptPat$alt)
imptPat$pos <- as.integer(imptPat$pos)

dndsskin_mySet = dndscv(imptPat, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
patsel_cv = dndsskin_mySet$sel_cv
patdNdS <- dndsskin_mySet$globaldnds

rm(substmodel,sel_cv,first_time,gr_genes,known_cancergenes,RefCDS)
#------Model Data------#
modelDat <- read.csv("~/Desktop/Darryl_collab/Model_Data_Analysis/HomeostaticDataAnalysis/OutputTables/Neutral/3D/vcfFilteredModelProbs/ModelAllPatients3DNeutral.txt",header=F,sep="\t")
colnames(modelDat) <- c("sampleID","chr","pos","ref","alt"  )
modelDat$sampleID <- as.character(modelDat$sampleID)
modelDat$chr <- as.character(modelDat$chr)
modelDat$ref <- as.character(modelDat$ref)
modelDat$alt <- as.character(modelDat$alt)
modelDat$pos <- as.integer(modelDat$pos)

dnds_model = dndscv(modelDat, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
modelsel_cv = dnds_model$sel_cv
modeldNdS <- dnds_model$globaldnds


#----Plot Results-----#
df <- read.csv("~/Desktop/dNdSModelToPatients.txt",header=F,sep="\t")
colnames(df) <- c("Source","Type","Type","mle","ci_low","ci_high")
df <- df[c(1,5,6,10),1:6]

ggplot(data=df, aes(x=Type,y=mle, fill=Source)) + geom_bar(position = "dodge", stat="identity") +
    theme_minimal(base_size = 14) + xlab("") + scale_fill_manual(values=c("darkgrey","lightgrey")) +
    ylab("dN / dS") + geom_errorbar(aes(ymin=ci_low, ymax=ci_high),
width=.1,                    # Width of the error bars
position=position_dodge(.9)) +
    geom_hline(mapping = NULL, data = NULL, yintercept=1, na.rm = FALSE,
    show.legend = NA, linetype = 2, color="red") +
    scale_y_continuous(breaks=c(0.0,0.5,1.0,1.5,2.0), limits=c(0.0,1.5))










