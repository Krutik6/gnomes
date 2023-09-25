#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#call libraries
library(httr2)
library(dplyr)
library(MultiAssayExperiment)
library(data.table)
library(gtools)
library(tidyr)
source("functions.R")

newdata <- read.csv("Trio7_forKrutik-nodups-060922.csv")
newdata <- newdata[-c(1)]
snvs <- newdata[which(newdata$VARIANT_CLASS == "SNV"),]
sep <- data.frame("Location" = snvs$Location)
sep <- as.data.frame(sep[!duplicated(sep$Location),])
sep <- separate(data = sep, col = 1, into = c("Chrom", "Pos"), sep = ":")
colnames(sep) <- c("Chrom", "Position")
sep$Type <- "SNV"
df <- sep
df$Chrom <- sub(df$Chrom, pattern = "chr", replacement = "")
rownames(df)<- NULL

df <- df[1:10,]
#add required nomenclature for further functions
df <- gnomeName(x = df, type = "SNV")
#make list of queries
queryMakerList <- gnomeQuest(x = df, type = "SNV")
#send list of queries to gnomad
querySentList <- gnomeCatch2(x = df, gnomeQuestList = queryMakerList,
                           batch = 30, sleep = 5, type = "SNV")
#check for errors
queryCheckList <- gnomeRatify(catchList = querySentList)
#if missing samples were found will need to re-run
missingList <- which(names(querySentList) %in% names(queryCheckList) == FALSE)
queryCheckList <- c(queryCheckList, querySentList[missingList])
#create list of variants
queryVariantsList <- gnomeVars(x = df, querySentList = queryCheckList,
                               type = "SNV")
#set up data for item retrieval
variants<- queryVariantsList[lapply(queryVariantsList,length)>0]
keepers <- names(variants)
siteData <- queryCheckList[keepers]
number <- length(keepers)

queryAlleleStats <- gnomeAlleles(variantList = variants,
                                 siteList = siteData)

gnomeDF <- gnomeDataFrame(listOfAlleles = queryAlleleStats, 
               variants = variants, siteData = siteData, 
               varNames = keepers, allVariants = queryVariantsList)
write.table(x = gnomeDF, "first800SNVs.csv")
remove(querySentList)
remove(queryCheckList)
################################################################################
newdata <- read.csv("Trio7_forKrutik-nodups-060922.csv")
newdata <- newdata[-c(1)]
indels <- newdata[which(newdata$Allele == "-"),]
sep <- data.frame(indels$Location)
colnames(sep) <- "Location"
sep <- separate(data = sep, col = Location, into = c("Chrom", "Pos"), sep = ":")
sep <- separate(data = sep, col = Pos, into = c("Start", "Stop"), sep = "-")
sep <- sep[!is.na(sep$Stop),]
sep$Type <- "INDEL"
df <- sep
df$Chrom <- sub(df$Chrom, pattern = "chr", replacement = "")

#test functions for indels
queryMakerList <- gnomeQuest(x = df, type = "INDEL")

querySentList <- gnomeCatch(x = df, gnomeQuestList = queryMakerList,
                            batch = 6, sleep = 20, type = "INDEL")

queryCheckList <- gnomeRatify(catchList = querySentList, type = "INDEL")

queryVariantsList <- gnomeVars(x = df, querySentList = queryCheckList,
                               type = "INDEL")

variants<- queryVariantsList[lapply(queryVariantsList,length)>0]
keepers <- names(variants)
siteData <- queryCheckList[keepers]
number <- length(keepers)

queryAlleleStats <- gnomeAlleles(variantList = variants,
                                 siteList = siteData)

gnomeDF <- gnomeDataFrame(listOfAlleles = queryAlleleStats, 
                          variants = variants, siteData = siteData, 
                          varNames = keepers, allVariants = queryVariantsList)

write.table(x = gnomeDF, "indels.csv")
