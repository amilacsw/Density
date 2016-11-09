d <- read.table("DataOut", header=TRUE)

library(ggplot2)

d$Variant <- factor(d$Variant, levels = d$Variant)

sp <- ggplot(d, aes(x=Variant,y=ClusterID)) + geom_point(aes(size=Probability))

sp+theme(axis.text.x=element_text(angle=90, size=8))

