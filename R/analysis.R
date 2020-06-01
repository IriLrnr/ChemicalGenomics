# include format_data
source("./R/format_data.R")

################# MAKE PRESENCE-ABSCENCE (OF ENRICHMENT) TABLE ###################
# Transform NA to 0 (just for it to be FALSE)
data.trans[is.na(data.trans)] <- 0
# Transform NA to 0 (just for it to be FALSE)
data.quanti[is.na(data.quanti)] <- 0
# save data
data.tf <- data.trans

### CUT 0
# Transform in presence-absence table (TRUE or FALSE)
data.tf <- data.tf[,] > 0
# Pra cada gene, quantos fármacos enriquecem
n.pergene <- apply(X = data.tf, MARGIN = 1, FUN = sum)
# Numero de genes enriquecidos por farmaci
n.permed <- apply(X = data.tf, MARGIN = 2, FUN = sum)
# Total of positive relations
total.enrichement <- sum(n.pergene)
##################################################################################


#################### TOP 10 ENRICHED TRANSPORTERS PER MED ###############################
# Create a table to store the sorted values
data.ordered.inc <- tibble(.rows = 121)
# Create columns orderd by enrichment, bind to table
for (i in 1:length(meds)) {
  order.inc <- as.factor(order(data.trans[,i], decreasing = TRUE))
  data.ordered.inc <- cbind(data.ordered.inc, order.inc)
}
# Rename columns to medsz
colnames(data.ordered.inc) <- meds
# Get the top10
data.top10 <- data.ordered.inc[1:10,]
# Create a new table to "translate" to gene names
data.t10.named <- tibble(.rows = 10)
# Fill the table
for (i in 1:length(meds)) {
  names <- paste(genes[data.top10[,i]])
  data.t10.named <- cbind(data.t10.named, names)  
}
# Rename columns to meds
colnames(data.t10.named) <- meds

data.t10.1<-data.t10.named[,1:500]
data.t10.2<-data.t10.named[,501:1000]
data.t10.3<-data.t10.named[,1001:1500]
data.t10.4<-data.t10.named[,1501:2000]
data.t10.5<-data.t10.named[,2001:2500]
data.t10.6<-data.t10.named[,2501:3000]
data.t10.7<-data.t10.named[,3001:3356]

write.csv(data.t10.named, "./output/top10_value.csv", row.names = FALSE, quote = FALSE)
##################################################################################

#################### BOTTOM 10 ENRICHED TRANSPORTERS PER MED ###############################
# Create a table to store the sorted values
data.ordered.dec <- tibble(.rows = 121)
# Create columns orderd by enrichment, bind to table
for (i in 1:length(meds)) {
  order.dec <- as.factor(order(data.trans[,i], decreasing = FALSE))
  data.ordered.dec <- cbind(data.ordered.dec, order.dec)
}
# Rename columns to medsz
colnames(data.ordered.dec) <- meds
# Get the top10
data.bottom10 <- data.ordered.dec[1:10,]
# Create a new table to "translate" to gene names
data.b10.named <- tibble(.rows = 10)
# Fill the table
for (i in 1:length(meds)) {
  names <- paste(genes[data.bottom10[,i]])
  data.b10.named <- cbind(data.b10.named, names)  
}
# Rename columns to meds
colnames(data.b10.named) <- meds

data.b10.1<-data.b10.named[,1:500]
data.b10.2<-data.b10.named[,501:1000]
data.b10.3<-data.b10.named[,1001:1500]
data.b10.4<-data.b10.named[,1501:2000]
data.b10.5<-data.b10.named[,2001:2500]
data.b10.6<-data.b10.named[,2501:3000]
data.b10.7<-data.b10.named[,3001:3356]

write.csv(data.b10.named, "./output/botton10_value.csv", row.names = FALSE, quote = FALSE)
##################################################################################

#################### UNIFY TOP AND BOTTON 10 #####################################
topbotton.1 <- rbind(data.t10.1, data.b10.1)
write.csv(topbotton.1, "./output/topbotton/topbotton_1.csv", row.names = FALSE, quote = FALSE)

topbotton.2 <- rbind(data.t10.2, data.b10.2)
write.csv(topbotton.2, "./output/topbotton/topbotton_2.csv", row.names = FALSE, quote = FALSE)

topbotton.3 <- rbind(data.t10.3, data.b10.3)
write.csv(topbotton.3, "./output/topbotton/topbotton_3.csv", row.names = FALSE, quote = FALSE)

topbotton.4 <- rbind(data.t10.4, data.b10.4)
write.csv(topbotton.4, "./output/topbotton/topbotton_4.csv", row.names = FALSE, quote = FALSE)

topbotton.5 <- rbind(data.t10.5, data.b10.5)
write.csv(topbotton.5, "./output/topbotton/topbotton_5.csv", row.names = FALSE, quote = FALSE)

topbotton.6 <- rbind(data.t10.6, data.b10.6)
write.csv(topbotton.6, "./output/topbotton/topbotton_6.csv", row.names = FALSE, quote = FALSE)

topbotton.7 <- rbind(data.t10.7, data.b10.7)
write.csv(topbotton.7, "./output/topbotton/topbotton_7.csv", row.names = FALSE, quote = FALSE)
##################################################################################
#################### TOP 10 ENRICHED GENES PER MED ###############################
# Create a table to store the sorted values
data.ordered <- tibble(.rows = 4810)
# Create columns orderd by enrichment, bind to table
for (i in 1:length(meds)) {
  order <- as.factor(order(data.quanti[,i], decreasing = TRUE))
  data.ordered <- cbind(data.ordered, order)
}
# Rename columns to medsz
colnames(data.ordered) <- meds
# Get the top10
data.top10 <- data.ordered[1:10,]
# Create a new table to "translate" to gene names
data.t10.named <- tibble(.rows = 10)
# Fill the table
for (i in 1:length(meds)) {
  names <- paste(genes[data.top10[,i]])
  data.t10.named <- cbind(data.t10.named, names)  
}
# Rename columns to meds
colnames(data.t10.named) <- meds

n.transporters <- vector()
#subset the transporters in the top10 table for the whole database
for (i in meds) {
  inters <- intersect(data.t10.named[,i], transps)
  n.transporters[i] <- length(inters)
}

# Count the number of meds that have transporters in the 10 most affected genes
transporters.tf <- n.transporters[] > 0
sum(transporters.tf)

write.csv(n.transporters, "./output/transporter_frequency.csv", row.names = TRUE, quote = FALSE)
##################################################################################

##################### 50 MOST FREQUENTLY ENRICHED GENES ##########################
# Create a table of occurences (in the top10 table) for each gene
occurences <- table(unlist(data.top10))
# Order the occurences to get the 50 genes that appear the most
ordered.occurence <- order(occurences, decreasing = TRUE)
# Sort occurences to get the frequencies
sorted.occurence <- sort(occurences, decreasing = TRUE)
ngenes <- sum(sorted.occurence > 0)
ngenes0 <- sum(sorted.occurence == 0)
# Extract the 50 first 
top50 <- ordered.occurence[1:50]
# Get the frequecies of the top 50 genes
top50.freq <- sorted.occurence[1:50]
# Take the names from the gene vector and put in the top 50
top50.names <- paste(meds[top50])
# Bind genes and frquency
top50.names.freq <- cbind(top50.names, top50.freq)
# Give names to columns
colnames(top50.names.freq) <- c("gene", "frequency")
# Write tables in csv
write.csv(top50.names.freq, "./output/freqtop50genes_trans.csv", row.names = FALSE, quote = FALSE)
####################################################################################

