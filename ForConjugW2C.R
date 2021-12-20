library (ggplot2)
library(readr)

data <- read_csv("C:/Users/marks/Desktop/Git/Project_Plasmid_C-/build/ForConjug2C3.csv")
                 
ggplot(data=data, aes(x= conjugationW)) +
  geom_point(aes(y = popsize, color = population, shape = location), size = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  
  labs(title = "The effect of increased conjugation at the wall \non equilibrium population sizes", 
       y = "population size (N0/N1)",
       x = "conjugation rate W")
