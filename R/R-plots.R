library(ggplot2)

d <- read.csv(file = "Model 0.1.csv", header = FALSE, sep = "\t")
d <- `colnames<-`(d, c("Generation", "P01", "Stdev1", "P02", "Stdev2", "P03", "Stdev3"))

ggplot(data = d, aes(x = Generation)) + 
  geom_point(aes(y = P01, colour = "High"), size = 2) + 
  geom_point(aes(y = P02, colour = "Low"), size = 2) + 
  geom_point(aes(y = P03, colour = "Medium"), size = 2) +
  labs(x = "Generation", y = "Mean P0", colour = "Cooperativeness") + 
  xlim(0, 10000) + ylim(0.0, 1.0)