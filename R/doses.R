# Load used libraries
library(openxlsx)

# Source analysis to make the tables
source("./R/analysis.R")

# Read coumpound information from HipHop library
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
transpose.t10 <- cbind(colnames(data.tb10.named), t(data.tb10.named[1:10,]))
colnames(transpose.t10)[1] <- c("Screen.ID")

compound.t10<- merge(compound.table[,1:4], transpose.t10, by="Screen.ID")
compound.t10 <- compound.t10[order(compound.t10$name, compound.t10$conc),]

# Print with colors to xlsx for Bessie
wb <- createWorkbook()
#add a worksheet to the workbook
addWorksheet(wb, "Compound_t10", gridLines = TRUE)
# write my analysis into the worksheet of the workbook, 
writeData(wb, "Compound_t10", compound.t10)
# here search for the respective HEX color-code and assign a name to the style
warm1Style <- createStyle(fontColour = "black", bgFill = "firebrick")

for(i in 1:nrow(transporters)){
  conditionalFormatting(wb, "Compound_t10", cols = 5:14,
                        rows = 1:nrow(compound.t10), rule = paste("=", transporters[i,1], sep = ""), style = warm1Style,
                        type = "contains")
}
saveWorkbook(wb, "./output/compounds_t10.xlsx", overwrite = TRUE)

# account the condition where "Significantly Warmer" is contained in a cell,
# then apply the respective style to it (in this case, warm1Style)
#############################################################################

