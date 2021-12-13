library (ggplot2)

data <- AlteringC

ggplot(data=data, aes(x=c)) +
  geom_line(aes(y = popsize, color = popid), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +

  labs(title = "The effect of conjugation on equilibrium population sizes", 
       y = "population size (N0/N1)",
       x = "conjugation rate")
