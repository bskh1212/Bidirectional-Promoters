##set working directory
setwd("~/Leeds/Project/5. GO term analysis")

##read list of temperature related GO terms
GO_temp <- read.delim("Temperature_GO_Terms.csv",
                      header=T,
                      sep=",")

##change column names
colnames(GO_temp) <- c("GO_Term", "Definition", "QuickGO_Search") 

##Get list of GO terms
temp_GO_terms <- GO_temp$GO_Term

##read list of bidirectional gene GO terms
BIP_GO <- read.delim("BIP_GO_terms.csv",
                     header=T,
                     sep=",")
##change column names
colnames(BIP_GO) <- c("Gene_ID", "GO_Term", "Definition",
                      "Accession", "Evidence_Code", "Domain")

##add column "Temp" saying if GO terms are temperature related or not 
BIP_GO$Temp <- BIP_GO$Accession %in% temp_GO_terms

##subset by temperature related GO terms
BIP_GO_temp <- BIP_GO[BIP_GO$Temp == T,]

##write temperature GO term genes to file
write.csv(BIP_GO_temp,
          "BIP_GO_temp.csv",
          row.names=F)

