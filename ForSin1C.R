library (ggplot2)
library(readr)

data <- results1comp0_35

ggplot(data=data, aes(x=pars)) +
  geom_line(aes(y = popsize, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  annotate("text", x=53, y=0.0001, label= "Flowrate = 0.35") +

  labs(title = "The effect of resource input concentration on\n equilibrium population sizes", 
       y = "Concentration of bacteria per ml",
       x = "resource input concentration (Sin)")
