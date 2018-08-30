library(tidyverse)

data.small <- data.table::fread("data/6eq_ssa_prob_extinction_small_R2.csv") %>%
	rename(gamma = V1, R1.init = V2, R2.init = V3, prob = V4) %>%
	filter(
		R2.init != 0,
		R1.init != 0
	) %>%
	group_by(gamma) %>%
	summarise(
		prob.mean = mean(prob),
		prob.sd = sd(prob)
	)

data.big <- data.table::fread("data/6eq_ssa_prob_extinction_fucking_big.csv") %>%
	rename(gamma = V1, R1.init = V2, R2.init = V3, prob = V4) %>%
	filter(
		R2.init != 0,
		R1.init != 0,
		R2.init >= 100
	) %>%
	group_by(gamma) %>%
	summarise(
		prob.mean = mean(prob),
		prob.sd = sd(prob)
	)

data.trans <- data.table::fread("data/6eq_ssa_prob_extinction_fucking_big.csv") %>%
	rename(gamma = V1, R1.init = V2, R2.init = V3, prob = V4) %>%
	filter(
		R2.init != 0,
		R1.init != 0,
		R2.init >= 30,
		R2.init <= 90
	) %>%
	group_by(gamma) %>%
	summarise(
		prob.mean = mean(prob),
		prob.sd = sd(prob)
	)

data.2eq <- data.table::fread("data/2eq_ssa_prob_extinction.csv") %>%
	# rename(gamma = V1, R1.init = V2, R2.init = V3, prob = V4) %>%
	# filter(
	# 	R2.init != 0,
	# 	R1.init != 0,
	# 	R2.init >= 30,
	# 	R2.init <= 90
	# ) %>%
	group_by(gamma) %>%
	summarise(
		prob.mean = mean(prob),
		prob.sd = sd(prob)
	)


plot_prob_exit_mean <- function(df) {
gg <- ggplot(df, aes(y = prob.mean, x = gamma)) +
	geom_errorbar(aes(ymin = prob.mean - prob.sd, ymax = prob.mean + prob.sd, color = "sd")) +
	geom_line(aes(color = "line"), linetype = "longdash", size = 0.5) +
	geom_point(aes(color = "mean")) +
	scale_x_continuous(breaks = seq(0, 0.55, 0.05)) +
	scale_y_continuous(limits = c(NA, NA), breaks = seq(0, 1, 0.1)) +
	labs(x = "Gamma", y = "Probability of Extinction") +
	scale_color_manual(name = NULL,
										 values = c("mean" = "steelblue1",
																 "sd" = "black",
																"line" = "grey35"
																),
										 breaks = c("mean", "sd"),
										 labels = c("Mean", "Standard deviation")) +
	# scale_linetype_manual(values = c("mean" = "solid"
	# 																 ),
	# 															guide = FALSE) +
	guides(color = guide_legend(override.aes = list(shape = c(19,NA),
																									linetype = c("blank", "solid")))) +
	theme_bw()

return(gg)
}

plot_prob_exit_mean(data.small) + theme_bw(base_family = "LM Roman 10") + alfR::ggexport(size = c(6, 2.5), "gg_6eq_plot_errorbars_small.pdf",  font = "lmodern")
plot_prob_exit_mean(data.trans) + theme_bw(base_family = "LM Roman 10") + alfR::ggexport(size = c(6, 2.5), "gg_6eq_plot_errorbars_trans.pdf",  font = "lmodern")
plot_prob_exit_mean(data.big) + theme_bw(base_family = "LM Roman 10") + alfR::ggexport(size = c(6, 2.5), "gg_6eq_plot_errorbars_big.pdf",  font = "lmodern")

plot_prob_exit_mean(data.2eq)
