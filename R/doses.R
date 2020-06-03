library(openxlsx)

source("./R/analysis.R")

compound.table <- read.csv("./data/raw/compound_table.csv", header = T, sep = ";")

###################### DIFFERENT DOSES ANALYSIS ############################
name.freq <- as.data.frame(table(compound.table$name))
name.freq <- subset(name.freq, Freq > 1)

inchikey.freq <- as.data.frame(table(compound.table$inchikey))
inchikey.freq <- subset(inchikey.freq, Freq > 1)

IUPAC.name.freq <- as.data.frame(table(compound.table$IUPAC.name))
IUPAC.name.freq <- subset(IUPAC.name.freq, Freq > 1)
IUPAC.name.freq <- IUPAC.name.freq[-1,]

doses <- subset(compound.table, compound.table$name %in% name.freq$Var1 &
                                compound.table$inchikey %in% inchikey.freq$Var1 &
                                compound.table$IUPAC.name %in% IUPAC.name.freq$Var1,
                select = c("Screen.ID", "name", "conc", "unit", "inchikey", "IUPAC.name"))

doses <- doses[order(doses$name, doses$inchikey, doses$IUPAC.name, doses$conc),]
doses.freq <- subset(as.data.frame(table(doses$name)), Freq > 1)

trans.rep <- subset(data.trans, select = c(paste(doses$Screen.ID)))
#############################################################################

####################### MERGE COMPOUND_TABLE AND TOP10 ######################
transpose.t10 <- cbind(rownames(transpose.t10), t(data.tb10.named[1:10,]))
colnames(transpose.t10)[1] <- c("Screen.ID")

compound.t10<- merge(compound.table[,1:4], transpose.t10, by="Screen.ID")
compound.t10 <- compound.t10[order(compound.t10$name, compound.t10$conc),]
#############################################################################

