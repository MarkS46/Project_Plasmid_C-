library (ggplot2)
library (ggpubr)

data <- conjsin2

data$pars2 <- as.factor(data$pars2)

fig_1 <- ggplot(data=data, aes(x=pars1, color = pars2)) +
  geom_line(aes(y = N0)) +
  scale_y_continuous(trans = "log10") +
  labs(title = "The equilibrium value of N0 for different conjugation rates\nand different input rates",
       y = "Population size",
       x = "Sin")

fig_2 <- ggplot(data=data, aes(x=pars1, color = pars2)) +
  geom_line(aes(y = N1)) +
  scale_y_continuous(trans = "log10") +
  labs(title = "The equilibrium value of N1 for different conjugation rates\nand different input rates",
       y = "Population size",
       x = "Sin")

ggarrange(fig_1, fig_2)

