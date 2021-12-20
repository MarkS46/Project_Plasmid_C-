library (ggplot2)
library(readr)

data <- read_csv("C:/Users/marks/Desktop/Git/Project_Plasmid_C-/build/ForSin1C8.csv")

ggplot(data=data, aes(x=Sin)) +
  geom_line(aes(y = popsize, color = population), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  annotate("text", x=53, y=0.0001, label= "Conjugation: 1e-8.8") +

  labs(title = "The effect of resource input on equilibrium population sizes", 
       y = "population size (N0/N1)",
       x = "input resource")
