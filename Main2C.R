library (ggplot2)
library(readr)

data <- Main2C

ggplot(data=data, aes(x= pars1)) +
  geom_point(aes(y = popsize, color = location, shape = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = location, linetype = population), size = 1, alpha =1)+
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") + 
  coord_cartesian(xlim = c(1e-6, NA)) +
  annotate("text", x=(1e-3 + 0.001 ), y=0.05, label= " cW = 1e-8.5" , size = 6) +

  labs(title = "The effect of attachment to the wall for N0\n and N1 on equilibrium population sizes (A2)", 
       y = "concentration of bacteria per ml",
       x = "attachment rate KLW1 (x 10 = KLW0)")+
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 10)
  )

