library (ggplot2)
library(readr)

data <- read_csv("C:/Users/marks/Desktop/Git/Project_Plasmid_C-/resultsKLW01.csv")

ggplot(data=data, aes(x= pars1)) +
  geom_point(aes(y = popsize, color = location, shape = population), size = 1.5, alpha = 1) +
  geom_line(aes(y = popsize, color = location, linetype = population), size = 1, alpha =1)+
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +  
  
  labs(title = "The effect of attachment to the wall for N0 & N1 \non equilibrium population sizes", 
       y = "population size (N0/N1)",
       x = "attachment rate KLW1 (KLW0 * 10)")
