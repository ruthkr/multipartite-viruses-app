library(tidyverse)

system("gcc -Wall -lm rk4_bipartite_4eq.c -o rk4-test")
system("./rk4-test 10 10 0 0 > results-test.csv")

data <- data.table::fread("results-test.csv")

ggplot(data) +
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "R2")) +
	geom_line(aes(x = V1, y = V4, colour = "p")) +
	geom_line(aes(x = V1, y = V5, colour = "S")) +
	labs(x = "Step", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R1" = "red3",        "R2" = "royalblue",
																"p" = "springgreen4", "S" = "purple3" )) +
	theme_bw()
