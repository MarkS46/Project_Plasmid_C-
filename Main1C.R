library (ggplot2)
library(readr)

data <- Main1AC9

ggplot(data=data, aes(x=pars)) +
  geom_point(aes(y = popsize, color = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  annotate("text", x=57, y=0.0001, label= "Ain = 0" , size = 8) +


  labs(title = "The effect of resource input concentration \non equilibrium population sizes", 
       y = "Concentration of bacteria per ml",
       x = "Resource input concentration (Sin)")  + 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13)
  )
