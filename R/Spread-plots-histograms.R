setwd("C:/Users/FSE-Education/Cooperation/Data/Model 1.1")

library(ggplot2)

filename <- "spread_model1.1_0.05.csv"

d <- read.csv(file = filename, header = FALSE, sep = "\t")

p <- ggplot(data = d, aes(x = d$V1))

for (i in 2:length(d[,1])) {
  write(paste0("p <- p + geom_point(aes(y = d$V", i, "))"), "tmp.R")
  source("tmp.R")
  file.remove("tmp.R")
}
p + ylim(0.0, 1.0) + ggtitle(filename)

d2 <- d[length(d[,1]),]
d2 <- as.data.frame(t(d2))
d2 <- as.data.frame(d2[-1,])
colnames(d2) <- "value"

ggplot(d2, aes(x = value)) + geom_histogram(bins = 100) + ggtitle(filename)
