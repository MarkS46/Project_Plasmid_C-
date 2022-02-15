library (ggplot2)
library(readr)

data <- Main1ACAdisp2

ggplot(data=data, aes(x=pars)) +
  geom_point(aes(y = popsize, color = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  annotate("text", x=2, y= 0.001, label= "Adisp1 = 1e-4" , size = 5) +


  labs(title = "The effect of Ain \non equilibrium population sizes", 
       y = "Concentration of bacteria per ml",
       x = "antibiotic input concentratoin (Ain)")  + 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13)
  )
