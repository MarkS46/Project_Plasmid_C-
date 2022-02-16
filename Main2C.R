library (ggplot2)
library(readr)

data <- Main2AC4

ggplot(data=data, aes(x= pars1)) +
  geom_point(aes(y = popsize, color = location, shape = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = location, linetype = population), size = 1, alpha =1)+
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") + 
  annotate("text", x=(1e-9 ), y=0.00005, label= "AinL = 10" , size = 6) +

  labs(title = "The effect of attachment to the wall for N1\non equilibrium population sizes", 
       y = "concentration of bacteria per ml",
       x = "attachment rate KLW1")+
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 10)
  )

