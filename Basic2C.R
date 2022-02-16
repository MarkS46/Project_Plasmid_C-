library (ggplot2)
library(readr)
library(ggpubr)

data <- nothing

data1 <- filter(data, location != "Wall")
data2 <- filter(data, location != "Lumen")

ggplot(data=data1, aes(x=time)) +
  geom_line(aes(y = size, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +

  
  labs(title = "The effect of plasmid on population \nconcentrations for 2 different compartments\n                      LUMEN",
       y = "concentration of bacteria in ml",
       x = "Time")+ 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13))


ggplot(data=data2, aes(x=time)) +
  geom_line(aes(y = size, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  
  
  labs(title = "The effect of plasmid on population \nconcentrations for 2 different compartments\n                      WALL",
       y = "concentration of bacteria in ml",
       x = "Time")+ 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13))
