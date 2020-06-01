## WORKING DIRECTORY SHOULD BE "Ludi compounds" ##
#library for plotting
library(ggplot2)

# Transform compound's files data into object ludi.compounds and ludi.compounds.sig 
source("./R/format_data_ludi.R")

# plot manhattan
# Adapted code for manhattan plot
axis.set <- ludi.compounds %>% 
  group_by(comp) %>% 
  summarize(center = (max(comp) + min(comp)) / 2)

ymin

ylim

ylim <- max(ludi.compounds$log2FoldChange) + 5
ymin <- min(ludi.compounds$log2FoldChange) - 5 
sig <- 1
sig2 <- -1
zero <- 0
nComp <- nrow(ludi.compounds)

breaks <- seq(-15, 15, 0.5)

labels <- as.character(breaks)
labels[!(breaks%%2==0)] <- ''
tick.sizes <- rep(.5, length(breaks))
tick.sizes[(breaks%%2==0)] <- 1

manhplot <- ggplot(ludi.compounds.nsig, aes(x = comp, y = log2FoldChange)) +
  geom_point(size = 1, alpha = 0.5, color = "darkslateblue") +
  geom_point(data = ludi.compounds.sig, size = 1, color = "firebrick3", alpha = 0.5, aes(x = comp, y = log2FoldChange)) +
  #geom_hline(yintercept = sig, color = "lightgoldenrod1", linetype = "dashed", size = 0.2) + 
  #geom_hline(yintercept = sig2, color = "lightgoldenrod1", linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = zero, color = "black", linetype = "solid", size = 0.2) + 
  scale_x_continuous(label = compounds, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), breaks = breaks, labels = labels, limits = c(ymin, ylim)) +
  ggtitle("Log2 FoldChange by compound") +
  labs(y = "Efluxers                          Influxers") + 
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust=1),
    axis.title.y = element_text(color="grey48", size=10),
    title = element_text(hjust = 0.5, color = "grey40", size = 10),
    plot.margin = margin(t = 3, r = 5, b = 3, l = 3, unit = "pt")
  )

manhplot

ggsave("./plots/manhattan_plot_certo.png", manhplot, height = 6, width = 5)

