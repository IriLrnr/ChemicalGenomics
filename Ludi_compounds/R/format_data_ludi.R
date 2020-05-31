library(tidyr)
library(dplyr)

# Set directory to compound file
setwd("./data/compound_files")
file.names <- dir()

# Create table named ludi.compounds, fill with info from the coumpound tables
ludi.compounds <- tibble()
remove.col <- c(9,8,5,4,2) 
length(file.names)
for (i in 1:length(file.names)) {
  temp <- read.table(file.names[i], header = T, sep = ",", quote = "\"", dec = ",")
  for (c in remove.col) {
    temp <- temp[,-c]  
  }
  comp <- rep(i, nrow(temp))
  temp <- cbind(temp, comp)
  ludi.compounds <- rbind(ludi.compounds, temp)
}

gsub(",", ".", ludi.compounds)

setwd("../../")

# Transform numbers in scientif form to doubles
ludi.compounds[,2] <-as.double(ludi.compounds[,2])

# Ignore lines with NA values for the numeric columns
ludi.compounds <- subset(ludi.compounds, log2FoldChange != is.na(log2FoldChange) &
                           pvalue != is.na(pvalue) &
                           padj != is.na(padj))

# Create a vector for significance
sig <- vector()
for (i in 1:nrow(ludi.compounds)) {
  if (ludi.compounds[i,3] <= 0.001 & ludi.compounds[i,4] <= 0.1) {
    sig[i] <- 1
  }
  else {
    sig[i] <- 0  
  }
}
# Bind the significance vector to the table of compounds
ludi.compounds <- cbind(ludi.compounds, sig)

# Create a second table with tha significant values only
ludi.compounds.sig <- subset(ludi.compounds, sig == 1)
ludi.compounds.nsig <- subset(ludi.compounds, sig == 0)

compounds <- read.csv("./data/Compostos")
compounds <- as.character(compounds[,1])

compounds[6] <- "Captan"
compounds[8] <- "Chlorothalonil"
compounds[11] <- "Difenoconazole"
compounds[13] <- "Epoxiconazole"
compounds[15] <- "Iprobenfos"
compounds[20] <- "Tebuconazole"
