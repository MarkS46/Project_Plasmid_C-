library (ggplot2)
library(readr)

data <- Main1ACH2

ggplot(data=data, aes(x=pars)) +
  geom_point(aes(y = popsize, color = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  annotate("text", x=5, y=1, label= "Ain = 20\nSin = 16\nd0/1 = 1e-2\nH0 = 4" , size = 5) +


  labs(title = "The effect of H1 \non equilibrium population sizes", 
       y = "Concentration of bacteria per ml",
       x = "half saturation constant H1 (plasmid bearing)")  + 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13)
  )
