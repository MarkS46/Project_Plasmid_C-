library (ggplot2)
library(readr)

data <- read_csv("C:/Users/marks/Desktop/Git/Project_Plasmid_C-/ForKLW02C2.csv")

ggplot(data=data, aes(x= AttachmentOfN0)) +
  geom_point(aes(y = popsize, color = location, shape = population), size = 1, alpha = 1) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  
  labs(title = "The effect of attachment at the wall for N1 \non equilibrium population sizes", 
       y = "population size (N0/N1)",
       x = "attachment rate KLW0")