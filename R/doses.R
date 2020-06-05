# Load used libraries
library(xlsx)

# Source analysis to make the tables
source("./R/analysis.R")

# Read coumpound information from HipHop library
compound.table <- read.csv("./data/raw/compound_table.csv", header = T, sep = ";")

###################### DOSES TABLE ############################
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
t10 <- as.data.frame(data.tb10.named[1:10,])

tf.t10 <- tibble(.rows = 10)
line <- vector()
for (i in 1:ncol(t10)) {
  line <- t10[,i] %in% transporters$V1
  tf.t10 <- cbind(tf.t10, line)
}

sum.t <-as.data.frame(as.factor(apply(tf.t10, MARGIN = 2, FUN = sum)))
colnames(sum.t) <- c("t")
transpose.t10 <- data.frame()
transpose.t10 <- cbind(colnames(t10), t(t10))
colnames(transpose.t10)[1] <- c("Screen.ID")

compound.t10<- merge(compound.table[,1:4], transpose.t10, by="Screen.ID")
compound.t10 <- compound.t10[order(compound.t10$name, compound.t10$conc),]
compound.t10 <- cbind(compound.t10, sum.t)

# Print with colors to xlsx for Bessie
cols <- length(compound.t10) - 1
sheetname <- "coumpound_t10"
write.xlsx(compound.t10, "./out_tables.xlsx", sheetName=sheetname, row.names = F)
file <- "out_tables.xlsx"
# but we want to highlight cells if value is equal to a transporter
wb <- loadWorkbook(file)              # load workbook
fo <- Fill(foregroundColor="lightgreen")   # create fill object
cs <- CellStyle(wb, fill=fo)        # create cell style
sheets <- getSheets(wb)               # get all sheets
sheet <- sheets[[sheetname]]          # get specific sheet
# get rows
rows <- getRows(sheet, rowIndex=1:nrow(compound.t10)+1)
# 1st row is headers
# get cells
cells <- getCells(rows, colIndex = 5:cols)
# extract the cell value
values <- lapply(cells, getCellValue)
highlight <- NULL
for (i in names(values)) {
  x <- as.factor(values[i])
  if (x %in% transporters$V1 && !is.na(x)) {
    highlight <- c(highlight, i)
  }
}
lapply(names(cells[highlight]), function(ii) setCellStyle(cells[[ii]], cs))
saveWorkbook(wb, file)
#####################################################################

####################### DOSES ANALYSIS ##############################
multiple.doses <- subset(compound.t10, compound.t10$name %in% doses$name, select = c("Screen.ID", "name", "conc", "unit", "t"))
sheetname = "multiple_doses"
write.xlsx(multiple.doses, "./out_tables.xlsx", sheetName=sheetname, row.names = F, append = TRUE)
#####################################################################
