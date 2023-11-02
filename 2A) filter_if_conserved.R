##Conserved Regions of Bidirectional Arrangement

##Take BIP regions from Arabidopsis and rice and see if they are present in
##wheat

##set working directory
setwd("~/Leeds/Project/3. Orthology")

##import orthologue data AT
##i.e. all of the arabidopsis bip genes that have wheat orthologues
orthologues <- read.delim("all_orthologues/AT_BIP_TA_orthologues.txt", sep=",",
                          header=T)
colnames(orthologues) <- c("Arabidopsis", "Triticum_aestivum", "Confidence")
orthologues <- na.omit(orthologues)

##import orthologue data OS
##i.e. all of the rice bip genes that have wheat orthologues
rice_orthologues <- read.delim("all_orthologues/OS_BIP_TA_orthologues.txt",
                               sep=",", header=T)
colnames(rice_orthologues) <- c("Oryza_sativa", "Triticum_aestivum",
                                "Confidence")
rice_orthologues <- na.omit(rice_orthologues)

##import putative BIP
BIP <- read.delim("../2. Head-to-Head Genes/BIPs/putative_BIP_triticum_aestivum.csv",
                  sep=",", header=T)

##keep the orthologues that are in the wheat BIP list only
orthologues <- orthologues[orthologues$Triticum_aestivum %in% BIP$Gene_ID,]
colnames(orthologues)[3] <- "AT_Confidence"

rice_orthologues <- rice_orthologues[rice_orthologues$Triticum_aestivum %in% BIP$Gene_ID,]
colnames(rice_orthologues)[3] <- "OS_Confidence"

##filter so there are only unique entries in the gene ID column. Sum the
##confidence values

##arabidopsis
repeated <- as.data.frame(table(orthologues$Triticum_aestivum))
repeated <- repeated[repeated$Freq>1,]$Var1
for (gene in repeated){
  value <- sum(orthologues[orthologues$Triticum_aestivum == gene, "AT_Confidence"])
  orthologues[orthologues$Triticum_aestivum == gene, "AT_Confidence"] <- value

}
orthologues <- unique(orthologues) 
colnames(orthologues)[3] <- "AT_Confidence_Sum"

length(unique(orthologues$Triticum_aestivum))

##rice
repeated <- as.data.frame(table(rice_orthologues$Triticum_aestivum))
repeated <- repeated[repeated$Freq>1,]$Var1
for (gene in repeated){
  value <- sum(rice_orthologues[rice_orthologues$Triticum_aestivum == gene, "OS_Confidence"])
  print(value)
  rice_orthologues[rice_orthologues$Triticum_aestivum == gene, "OS_Confidence"] <- value
  
}
rice_orthologues <- unique(rice_orthologues)
colnames(rice_orthologues)[3] <- "OS_Confidence_Sum"

##crossover between rice/arabidopsis/wheat
##how many BIP genes are conserved across all three? 154
table(unique(orthologues$Triticum_aestivum) %in% unique(rice_orthologues$Triticum_aestivum))
table(unique(rice_orthologues$Triticum_aestivum) %in% unique(orthologues$Triticum_aestivum))

##how many times does each BIP orthologue in triticum aestivum
##occur as an orthologue of an arabidopsis BIP gene?
ortho_number <- as.data.frame(table(orthologues$Triticum_aestivum))
colnames(ortho_number)[2] <- "AT_Freq"
table(ortho_number$AT_Freq)

##of a rice BIP gene?
rice_ortho_number <- as.data.frame(table(rice_orthologues$Triticum_aestivum))
colnames(rice_ortho_number)[2] <- "OS_Freq" 
table(rice_ortho_number$OS_Freq)

##merge BIP dataframe with AT orthologue information
BIP_ortho <- merge(BIP, orthologues, by.x="Gene_ID", by.y="Triticum_aestivum", all.x=T)
BIP_ortho$Arabidopsis <- NULL
BIP_ortho <- merge(BIP_ortho, ortho_number, by.x="Gene_ID", by.y="Var1", all.x=T)

##with OS orthologue info
BIP_ortho <- merge(BIP_ortho, rice_orthologues, by.x="Gene_ID", by.y="Triticum_aestivum", all.x=T)
BIP_ortho$Oryza_sativa <- NULL
BIP_ortho <- merge(BIP_ortho, rice_ortho_number, by.x="Gene_ID", by.y="Var1", all.x=T)

##remove duplicate wheat genes
BIP_ortho <- unique(BIP_ortho)

##start filtering by orthologue information

##define a function that removes unpaired genes from a dataset
##data frame must have column "PairID". returns dataframe
filter_unpaired <- function(df){
  pairs <- as.data.frame(table(df$PairID))
  pairs <- pairs[pairs$Freq==2,]$Var1
  df <- df[df$PairID %in% pairs,]
  return(df)
}

################################################################################

subset1 <- BIP_ortho

##How many genes are conserved as pairs? 44
subset2 <- subset1[subset1$AT_Freq>0,]
subset2 <- subset2[subset2$OS_Freq>0,]
subset2 <- na.omit(subset2)
subset2 <- filter_unpaired(subset2)

##How many genes are conserved as pairs in Arabidopsis? 106
subset2 <- BIP_ortho[!is.na(BIP_ortho$AT_Confidence_Sum),]
subset2 <- filter_unpaired(subset2)

##How many genes are conserved as pairs in rice? 518
subset3 <- BIP_ortho[!is.na(BIP_ortho$OS_Confidence_Sum),]
subset3 <- filter_unpaired(subset3)

##How many pairs have orthologues in rice AND Arabidopsis?
subset4 <- BIP_ortho[!is.na(BIP_ortho$AT_Confidence_Sum),]
subset4 <- subset4[!is.na(subset4$OS_Confidence_Sum),]
subset4 <- filter_unpaired(subset4)


##22 conserved pairs between wheat, rice, and arabidopsis
##none for either if you only count high confidence orthologues


