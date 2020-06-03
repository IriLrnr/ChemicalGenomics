source("./R/analysis.R")

compound.table <- read.csv("./data/raw/compound_table.csv", header = T, sep = ";")

name.freq <- as.data.frame(table(compound.table$name))
name.freq <- subset(name.freq, Freq > 1)

inchikey.freq <- as.data.frame(table(compound.table$inchikey))
inchikey.freq <- subset(inchikey.freq, Freq > 1)

IUPAC.name.freq <- as.data.frame(table(compound.table$IUPAC.name))
IUPAC.name.freq <- subset(IUPAC.name.freq, Freq > 1)
IUPAC.name.freq <- IUPAC.name.freq[-1,]

PCID.freq <- as.data.frame(table(compound.table$PCID))
PCID.freq <- subset(PCID.freq, Freq > 1)
PCID.freq <- PCID.freq[-1,]


doses <- subset(compound.table, compound.table$name %in% name.freq$Var1 |
                                compound.table$inchikey %in% inchikey.freq$Var1 |
                                compound.table$IUPAC.name %in% IUPAC.name.freq$Var1 |
                                compound.table$PCID %in% PCID.freq$Var1, 
                select = c("Screen.ID", "name", "conc", "unit", "PCID", "inchikey", "IUPAC.name"))

doses <- doses[order(doses$name, doses$inchikey, doses$PCID, doses$IUPAC.name),]
