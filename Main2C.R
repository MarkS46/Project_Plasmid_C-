library (ggplot2)
library(readr)

data <- fuckoffplz8

ggplot(data=data, aes(x= pars1)) +
  geom_point(aes(y = popsize, color = location, shape = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = location, linetype = population), size = 1, alpha =1)+
  scale_x_continuous(trans= "log10") +

 

  labs(title = "The effect of migratino and \nantibiotics on equilibrium population sizes", 
       y = "concentration of bacteria per ml",
       x = "KLW1 migration rate")+
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 10)
  )

