#Clear out leftover stuff from last session
rm(list=ls())

#########################################################################################
#install dependencies

#bioconductor and bioconductor packages
if (! "BiocInstaller" %in% rownames(installed.packages())){
  #Check if bioconductor is already installed
  source("https://bioconductor.org/biocLite.R")
  biocLite()
}

bioconductor_packages<-c("GenomicFeatures", "AnnotationDbi","BSgenome","BSgenome.Hsapiens.UCSC.hg19")
new.BC.packages<-bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
if (length(new.BC.packages)) {biocLite(new.BC.packages)}
#Normal R packages
list_of_packages<-c("vegan","RColorBrewer","gplots","ggplot2","data.table","parallel")
new.packages<-list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) {install.packages(new.packages)}

#Load dependencies

for (package in list_of_packages){library(package,character.only = TRUE)}

for (package in bioconductor_packages){library(package,character.only = TRUE)}

####################################################################################


#Setup variables for current working machine

home_dir=dirname(path.expand("~"))

OS_type=.Platform$OS.type

if (OS_type=="windows"){
  g_drive="Google Drive"
} else if (OS_type=="unix"){
  g_drive="gdrive"
}

project_name="20162017SumSemResearch"

project_dir=file.path(home_dir,g_drive,project_name)

data_directory=file.path(project_dir,"Datasets")

#Loading in data
SNP_data_dir=file.path(data_directory,"SNPDatabases")

GWAS_cat_SNPs=file.path(SNP_data_dir,"gwas_catalog_v1.0-associations_e86_r2016-10-23.tsv")

SNP_data=fread(GWAS_cat_SNPs,sep="\t")

SNP_data=clean_SNP_data(SNP_data)
#Imputing sequences for all of the indiviudals
#I Think I've misunderstood and this is not necessary

#Extracting flanking sequences
flanking_sequences=list()
flanking_offset=100

for (i in 1:nrow(SNP_data)) {
  SNP=SNP_data[i,]
  chromosomes=trimws(strsplit(SNP$CHR_ID,";")[[1]])
  positions=trimws(strsplit(SNP$CHR_POS,";")[[1]])
  SNP_names=trimws(strsplit(SNP$SNPS,";")[[1]])
  for (j in 1:length(chromosomes)){
    chr=paste("chr",chromosomes[j],sep="")
    pos=positions[j]
    SNP_name=SNP_names[j]
    flanking_sequences[[SNP_name]]=paste(getSeq(Hsapiens,chr,position-flanking_offset,position),
                                         getSeq(Hsapiens,chr,position+1,position+flanking_offset),
                                         sep='')
  }
}

#Example linked to by Loic
chr <- 'chr19'
position <- 59900243
alleles <- '[T/C]'
offset <- 60

test_sequence <- paste(getSeq(Hsapiens,chr,position-offset,position-1),
             alleles,
             getSeq(Hsapiens,chr,position+1,position+offset),
             sep='')
#Motif discovery on flanking sequences for a specific SNP

#Writing of intermediate file

#Collation of results across entire genome for background dataset

#Plotting of relative frequency of each motif around a SNP

#Plotting of enrichment of motif compared to entire background dataset


#########################################################################
#Functions
#########################################################################

clean_SNP_data<-function(SNP_data,exclude_interations){
  '''A function which takes in the original SNP data and then simplify the variable
  SNP naming convenients to always have a chromosome ID and position
  It also flattens data.table while splitting based on the ; character in the chromosome,positions
  and SNP_name columns.
  
  Input: SNP_data - data.table
         A data.table of the SNPS from the GWAS catalogue
  Output: cleaned_SNP_data - data.table
         A flattened and less ambiguous form of the GWAS catalogue'''
  SNP_data$UNIQ_ID=(1:nrow(SNP_data))
  if (exclude_interactions){
    working_SNP_data=SNP_data[!(grepl("x",SNP_data$CHR_ID) | grepl("x",SNP_data$CHR_POS)),]
  }
  else {
    working_SNP_data=SNP_data
  }
  #Splitting of SNPs of form chrID:Positions
  split_data=working_SNP_data[grepl("chr[0-9]{1,}:[0-9]{1,}",working_SNP_data$SNPS),]
  #Needs confirmation as to whether the extra characters here are intential or a mistake
  split_data$SNPS=gsub("[:-][A-Z]$","",split_data$SNPS)
  #I'll just exclude these if I can't confirm there properties.
  
  
  
  
  
  #splitting SNPs of form
  

  #Extract the problem sequences with no values
  missing_info=SNP_data[SNP_data$CHR_ID=="" | SNP_data$CHR_POS=="",]
}


out <- working_SNP_data[, list(CHR_ID=unlist(strsplit(CHR_ID, "(;|; )")),CHR_POS=unlist(strsplit(CHR_POS, "(;|; )")),SNPS = unlist(strsplit(SNPS, "(;|; )")),"STRONGEST SNP-RISK ALLELE" = unlist(strsplit(`STRONGEST SNP-RISK ALLELE`, "(;|; )"))), by=UNIQ_ID]
test_out<-working_SNP_data[,list(CHR_ID=lengths(regmatches(CHR_ID, gregexpr("(;|; )", CHR_ID))),CHR_POS=lengths(regmatches(CHR_POS, gregexpr("(;|; )", CHR_POS))),"STRONGEST SNP-RISK ALLELE"=lengths(regmatches(`STRONGEST SNP-RISK ALLELE`, gregexpr("(;|; )", `STRONGEST SNP-RISK ALLELE`))),SNPS=lengths(regmatches(SNPS, gregexpr("(;|; )", SNPS)))),by=UNIQ_ID]

out_new<-working_SNP_data[out$UNIQ_ID,]



#First
#Imputing sequences for all of the indiviudals
#I Think I've misunderstood and this is not necessary

#Extracting flanking sequences

#Example linked to by Loic
chr <- 'chr19'
position <- 59900243
alleles <- '[T/C]'
offset <- 60

test_sequence <- paste(getSeq(Hsapiens,chr,position-offset,position-1),
                       alleles,
                       getSeq(Hsapiens,chr,position+1,position+offset),
                       sep='')
#Motif discovery on flanking sequences for a specific SNP

#Writing of intermediate file

#Collation of results across entire genome for background dataset

#Plotting of relative frequency of each motif around a SNP

#Plotting of enrichment of motif compared to entire background dataset
