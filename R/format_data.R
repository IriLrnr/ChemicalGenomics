# Library
library(tidyr)
library(dplyr)

# Read HOP table
data.raw <- read.table("./data/raw/fitness_defect_matrix_hom.txt", header = TRUE, sep="\t", nrows=4811)
# Read transporters names:
transporters <-read.csv("./data/raw/transporters.csv", header = FALSE, sep= ",", nrows=130)

# Save raw
data.quanti <- data.raw


# Have to change name of the transporters column
data.trans <- filter(data.raw, (X %in% transporters$V1))

# Transform first row in line names in data.trans
rownames(data.trans) <- data.trans[,1]
data.trans <- subset(data.trans, select = -X)
# Transform first row in line names in data.quanti
rownames(data.quanti) <- data.quanti[,1]
data.quanti <- subset(data.quanti, select = -X)
# Make a vector with rownames and colnames
meds <- colnames(data.trans)
transps <- rownames(data.trans)
transp.hiphop <- rownames(data.trans);
genes <- rownames(data.quanti)

not.in <- subset(transporters, !(V1 %in% transp.hiphop))
write.csv(not.in, "./output/not_in_hop.csv", row.names = F, quote = F)
