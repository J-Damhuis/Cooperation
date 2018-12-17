library(ggplot2)

d <- read.csv(file = "R/Model 1.1.csv", header = FALSE, sep = "\t")

ggplot(data = d, aes(x = d[1])) + geom_point(aes(y = d[2], colour = "P0 = 0.95"), size = 0.5) + 
  geom_point(aes(y = d[6], colour = "P0 = 0.05"), size = 0.5) + 
  geom_point(aes(y = d[10], colour = "P0 = 0.67"), size = 0.5) +
  labs(x = "Generation", y = "Fraction of interactions where focal individual cooperates") +
  scale_colour_manual(values = c("orange", "grey", "blue")) +
  xlim(0, 1000000) + ylim(0.0, 1.0) + theme(legend.title=element_blank())

ggplot(data = d, aes(x = d[1])) + geom_point(aes(y = d[3], colour = "P0 = 0.95"), size = 0.5) + 
  geom_point(aes(y = d[7], colour = "P0 = 0.05"), size = 0.5) + 
  geom_point(aes(y = d[11], colour = "P0 = 0.67"), size = 0.5) +
  labs(x = "Generation", y = "Mean P0 value") +
  scale_colour_manual(values = c("orange", "grey", "blue")) +
  xlim(0, 1000000) + ylim(0.0, 1.0) + theme(legend.title=element_blank())

ggplot(data = d, aes(x = d[1])) + geom_point(aes(y = d[5], colour = "P0 = 0.95"), size = 0.5) + 
  geom_point(aes(y = d[9], colour = "P0 = 0.05"), size = 0.5) + 
  geom_point(aes(y = d[13], colour = "P0 = 0.67"), size = 0.5) +
  labs(x = "Generation", y = "Fraction of population which obtains information") +
  scale_colour_manual(values = c("orange", "grey", "blue")) +
  xlim(0, 1000000) + ylim(0.0, 1.0) + theme(legend.title=element_blank())
