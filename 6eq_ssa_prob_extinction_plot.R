library(tidyverse)

# source("6eq_find_first_point_to_P3.R")

results.df <- rbind(
	data.table::fread("data/6eq_ssa_prob_extinction_fucking_big.csv"),
	data.table::fread("data/6eq_ssa_prob_extinction_small_R2.csv"))%>%
	rename(gamma = V1, R1.init = V2, R2.init = V3, prob = V4)

results.df %>% filter(
	gamma == 0.52,
	R2.init == 20
	)

plot_ssa_prob_extinction <- function(df, seq, facet.wrap = F) {
	gg <- ggplot(df %>% filter(R2.init != 0, R1.init != 0, R2.init %in% seq)) +
		aes(x = gamma, y = prob, color = as.factor(R2.init)) +
		labs(color = "R2.init") +
		scale_x_continuous(breaks = seq(0, 0.55, 0.05)) +
		scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, 0.1)) +
		labs(x = "Gamma", y = "Probability of extinction", color = "Initial R2") +
		theme_bw()

	if(facet.wrap) {
		gg <- gg +
			geom_point(size = 0.5) +
			geom_line(linetype = "dashed", size = 0.5) +
			facet_wrap(~R2.init, ncol = 4) +
			scale_x_continuous(breaks = seq(0, 0.55, 0.2)) +
			scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, 0.2))
	} else {
		gg <- gg +
			geom_point() +
			geom_line(linetype = "dashed")
	}

	return(gg)
}

# New plots
plot_ssa_prob_extinction(results.df, seq(1,20,1))
plot_ssa_prob_extinction(results.df, seq(1,20,1), facet.wrap = T)


# Old plots
plot_ssa_prob_extinction(results.df, seq(0,100,1)) #+ theme_bw(base_family = "LM Roman 10") + alfR::ggexport("6eq_prob_exit_small.pdf", font = "lmodern")

plot_ssa_prob_extinction(results.df, seq(0,1000,100)) #+ theme_bw(base_family = "LM Roman 10") + alfR::ggexport("6eq_prob_exit_big.pdf", font = "lmodern")

plot_ssa_prob_extinction(results.df, seq(0,1000,50), facet.wrap = T) #+ theme_bw(base_family = "LM Roman 10") + theme(strip.background = element_blank(), strip.text.x = element_blank())+ alfR::ggexport(size = c(6,6), "6eq_prob_exit_big_facet.pdf", font = "lmodern")
