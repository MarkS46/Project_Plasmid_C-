library (ggplot2)
library(readr)
library(ggpubr)

data <- Basic2C

ggplot(data=data, aes(x=t)) +
  geom_line(aes(y = popsize, color = location, linetype = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  labs(title = "The effect of plasmid on population \nconcentrations for 2 different compartments",
       y = "concentration of bacteria in ml",
       x = "Time")+ 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13))



