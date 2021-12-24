library (ggplot2)
library(readr)
data <- Basic1C

data$`resource input` <- factor(data$`resource input`)

ggplot(data=data, aes(x=t)) +
  geom_line(aes(y = popsize, color = `resource input`, linetype = population), size = 1) +
  scale_color_manual(breaks = c("16"),
                    values=c("red")) +
  scale_y_continuous(trans = "log10") +
  labs(title = "The effect of plasmid on population\nconcentrations for resource input of 16",
       y = "concentration of bacteria in ml",
       x = "Time") + 
  
  theme(
    plot.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    axis.text = element_text(size= 13))

 



