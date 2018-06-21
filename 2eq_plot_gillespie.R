library(tidyverse)

data.ssa <- data.table::fread("data/results-ssa.csv")

ggplot(data.ssa) +
	geom_line(aes(x = V1, y = V2, colour = "R")) +
	geom_line(aes(x = V1, y = V3, colour = "p")) +
	labs(x = "Time", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R" = "red3",
																"p" = "springgreen4")) +
	theme_bw()
