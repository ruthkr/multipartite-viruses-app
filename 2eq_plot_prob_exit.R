library(tidyverse)

# Read data
data.first <- as.tibble(data.table::fread("data/ssa-prob-exit-results_final.csv"))
data.second <- as.tibble(data.table::fread("data/ssa-prob-exit-results.csv"))

# Compute the mean value and standard deviation
data.first <- data.first %>%
	dplyr::rename(gamma = V1, prob = V2) %>%
	group_by(gamma) %>%
	summarise(
		mean.prob = mean(prob),
		sd.prob = sd(prob)
	)

data.second <- data.second %>%
	dplyr::rename(gamma = V1, prob = V2) %>%
	group_by(gamma) %>%
	summarise(
		mean.prob = mean(prob),
		sd.prob = sd(prob)
	)

# Plot the probability for each gamma
ggplot(data.first) +
	aes(x = gamma, y = mean.prob) +
	geom_errorbar(aes(ymin = mean.prob - sd.prob, ymax = mean.prob + sd.prob, colour = "First")) +
	geom_point(aes(colour = "First")) +
	geom_line(linetype = "dashed", aes(colour = "First")) +
	geom_errorbar(data = data.second, aes(ymin = mean.prob - sd.prob, ymax = mean.prob + sd.prob, colour = "Second")) +
	geom_point(data = data.second, aes(colour = "Second")) +
	geom_line(data = data.second, linetype = "dashed", aes(colour = "Second")) +
	# stat_smooth(method = "loess") +
	scale_color_manual(values = c("First" = "maroon3",
																"Second" = "steelblue2"))+
	scale_x_continuous(breaks=seq(0,0.55,0.05))+
	labs(y = "Probability of Extinction", x="Value of Gamma", color = "Simulation")+
	theme_bw()

