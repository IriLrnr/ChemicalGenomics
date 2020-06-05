# include format_data
source("./R/format_data.R")

################# MAKE PRESENCE-ABSCENCE (OF ENRICHMENT) TABLE ###################
# Transform NA to 0 (just for it to be FALSE)
data.trans[is.na(data.trans)] <- 0
# Transform NA to 0 (just for it to be FALSE)
data.quanti[is.na(data.quanti)] <- 0
# save data
data.tf <- data.trans
##################################################################################

############### TOP AND BOTTOM 10 ENRICHED TRANSPORTERS PER MED ##################
# Create a table to store the sorted values
data.ordered <- tibble(.rows = 121)
# Create columns orderd by enrichment, bind to table
for (i in 1:length(meds)) {
  order <- as.factor(order(data.trans[,i], decreasing = FALSE))
  data.ordered <- cbind(data.ordered, order)
}
# Rename columns to medsz
colnames(data.ordered) <- meds
# Get the top10
trans.tb10 <- rbind(data.ordered[1:10,], data.ordered[112:121,])
# Create a new table to "translate" to gene names
trans.tb10.named <- tibble(.rows = 20)
# Fill the table
for (i in 1:length(meds)) {
  names <- paste(genes[trans.tb10[,i]])
  trans.tb10.named <- cbind(trans.tb10.named, names)  
}
# Rename columns to meds
colnames(trans.tb10.named) <- meds
##################################################################################

############## TOP AND BOTTOM 10 ENRICHED GENES PER MED ##########################
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
data.tb10 <- rbind(data.ordered[1:10,], data.ordered[3347:3356,])
# Create a new table to "translate" to gene names
data.tb10.named <- tibble(.rows = 20)
# Fill the table
for (i in 1:length(meds)) {
  names <- paste(genes[data.tb10[,i]])
  data.tb10.named <- cbind(data.tb10.named, names)  
}
# Rename columns to meds
colnames(data.tb10.named) <- meds

topbottom.1<-data.tb10.named[,1:500]
topbottom.2<-data.tb10.named[,501:1000]
topbottom.3<-data.tb10.named[,1001:1500]
topbottom.4<-data.tb10.named[,1501:2000]
topbottom.5<-data.tb10.named[,2001:2500]
topbottom.6<-data.tb10.named[,2501:3000]
topbottom.7<-data.tb10.named[,3001:3356]

write.csv(topbottom.1, "./output/topbottom/topbottom_1.csv", row.names = FALSE, quote = FALSE)
write.csv(topbottom.2, "./output/topbottom/topbottom_2.csv", row.names = FALSE, quote = FALSE)
write.csv(topbottom.3, "./output/topbottom/topbottom_3.csv", row.names = FALSE, quote = FALSE)
write.csv(topbottom.4, "./output/topbottom/topbottom_4.csv", row.names = FALSE, quote = FALSE)
write.csv(topbottom.5, "./output/topbottom/topbottom_5.csv", row.names = FALSE, quote = FALSE)
write.csv(topbottom.6, "./output/topbottom/topbottom_6.csv", row.names = FALSE, quote = FALSE)
write.csv(topbottom.7, "./output/topbottom/topbottom_7.csv", row.names = FALSE, quote = FALSE)

n.transporters <- vector()
#subset the transporters in the top10 table for the whole database
for (i in meds) {
  inters <- intersect(data.tb10.named[,i], transps)
  n.transporters[i] <- length(inters)
}

# Count the number of meds that have transporters in the 10 most affected genes
transporters.tf <- n.transporters[] > 0
sum(transporters.tf)

write.csv(n.transporters, "./output/transporter_frequency.csv", row.names = TRUE, quote = FALSE)
##################################################################################

##################### 50 MOST FREQUENTLY ENRICHED GENES ##########################
# Create a table of occurences (in the top10 table) for each gene
occurences <- table(unlist(data.tb10))
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

############################### RATIOS ##########################################
tb10.freq <- as.data.frame(table(unlist(data.tb10.named)))
transps.freq.10 <- subset(tb10.freq, Var1 %in% transps)
sum (transps.freq.10$Freq) / sum (tb10.freq$Freq)

tx.freq <- as.data.frame(table(unlist(data.tb10.named[1:1,])))
transps.freq.tx <- subset(tx.freq, Var1 %in% transps)
sum (transps.freq.tx$Freq) / sum(tx.freq$Freq)

bx.freq <- as.data.frame(table(unlist(data.tb10.named[29:20,])))
transps.freq.bx <- subset(bx.freq, Var1 %in% transps)
sum (transps.freq.bx$Freq) / sum(bx.freq$Freq)
###################################################################################
