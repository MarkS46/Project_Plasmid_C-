library (ggplot2)
library(readr)
data <- read_csv("C:/Users/marks/Desktop/Git/Project_Plasmid_C-/Basic1C.csv")

data$`resource input` <- factor(data$`resource input`)
data$trans <- log10(data$popsize)

ggplot(data=data, aes(x=t)) +
  geom_line(aes(y = popsize, color = `resource input`, linetype = population), size = 1) +
  scale_color_manual(breaks = c("50"),
                    values=c("red")) +
  scale_y_continuous(trans = "log10") +
  labs(title = "The effect of plasmid on growth rate \nfor resource input of 50",
       y = "Population size",
       x = "Time") 

 



