library(ggplot2)

source("./R/analysis.R")

limit <- length(meds)
 
beg <- 1 

# Create table for manhattan plot
man.compounds <- tibble(.rows = 16142360)
temp <- tibble()
strains <- 1:length(genes)
for (m in beg:limit){
  compound <- rep(m, length(strains))
  temp <- cbind(strains, compound, data.quanti[,m])
  man.compounds <- rbind(man.compounds, temp)
} #slow

# Eliminate NA and 0 rows
#man.compounds[is.na(man.compounds)] <- 0
#man.compounds.pos <- subset(man.compounds, V3 < 0, select = c("strains", "compound", "V3"))

man.compounds.pos <- subset(man.compounds, V3 > 5, select = c("strains", "compound", "V3"))
man.compounds.neg5 <- subset(man.compounds, V3 < -5, select = c("strains", "compound", "V3"))
man.compounds.pos <- rbind(man.compounds.pos, man.compounds.neg5)

# Make a second table for transporters
man.transport <- tibble()
temp <- tibble()
strains <- 1:length(transps)
for (m in beg:limit){
  compound <- rep(m, length(strains))
  temp <- cbind(strains, compound, data.trans[,m])
  man.transport <- rbind(man.transport, temp)
}
# Eliminate negative rows
man.transport[is.na(man.transport)] <- 0
#man.transport.pos <- subset(man.transport, V3 < 0, select = c("strains", "compound", "V3"))

man.transport.pos <- subset(man.transport, V3 > 5, select = c("strains", "compound", "V3"))
man.transport.neg5 <- subset(man.transport, V3 < -5, select = c("strains", "compound", "V3"))
man.transport.pos <- rbind(man.transport.pos, man.transport.neg5)

# Adapted code for manhattan plot
ylim <- max(man.compounds.pos$V3)
ymin <- min(man.compounds.pos$V3)
sig <- 5
sig2 <- -5
zero <- 0
nComp <- nrow(man.compounds.pos)

manhplot <- ggplot(man.compounds.pos, aes(x = compound, y = V3)) +
  geom_point(size = 0.3, alpha = 0.75, color = "lightseagreen") +
  geom_point(data = man.transport.pos, size = 0.3, color = "tomato", alpha = 1, aes(x = compound, y = V3)) +
  geom_hline(yintercept = sig, color = "lightgoldenrod1", linetype = "dashed") + 
  #geom_hline(yintercept = sig2, color = "lightgoldenrod1", linetype = "dashed") + 
  geom_hline(yintercept = zero, color = "grey40", linetype = "solid") + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(ymin, ylim)) +
  scale_size_continuous(range = c(0.5,3)) +
  #ggtitle("") +
  labs(x = "Compound", 
       y = "Fitness") + 
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = 'black'),
    axis.title.x = element_text(color="grey40", size=10),
    axis.text.x = element_blank(),
    axis.title.y = element_text(color="grey40", size=10),
    #title = element_text(hjust = 0.5, color = "grey40", size = 10),
    plot.margin = margin(t = 3, r = 5, b = 3, l = 3, unit = "pt")
  )

manhplot

ggsave("manhattan_maismenos5.png", manhplot)

