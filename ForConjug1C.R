library (ggplot2)
library(readr)

data <- read_csv("C:/Users/marks/Desktop/Git/Project_Plasmid_C-/ForConjug1C.csv")

ggplot(data=data, aes(x=c)) +
  geom_line(aes(y = popsize, color = popid), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +

  labs(title = "The effect of conjugation on equilibrium population sizes", 
       y = "population size (N0/N1)",
       x = "conjugation rate")
