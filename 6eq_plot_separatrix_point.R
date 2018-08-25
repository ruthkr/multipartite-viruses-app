library(tidyverse)

source("6eq_find_first_point_to_P3.R")

scale <- 1000

df.gg.sep <- df.clean.sep %>%
	filter(gamma %in% c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53)) %>%
	# filter(gamma %in% c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53)) %>%
	mutate(R1.init = R1.init * scale,
				 R2.init = R2.init * scale) %>%
	group_by(gamma, R1.init) %>%
	summarise(R2.init = min(R2.init))


df.gg.init <- df.clean.init %>%
	filter(gamma %in% c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53)) %>%
	# filter(gamma %in% c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53)) %>%
	mutate(R1.init = R1.init * scale,
				 R2.init = R2.init * scale)

plot_6eq_separatrix <- function(df.init, df.sep, drift = c(0, 0), limits = c(NA, NA), zoom.box = F) {
	gg <- ggplot(df.init) +
		# Initial condition
		geom_point(aes(x = (R1.init + drift[1]) , y = R2.init, color = as.factor(gamma))) +
		# Separatrix
		geom_line(data = df.sep, aes(x = (R1.init + drift[2]), y = R2.init, colour = as.factor(gamma)), linetype = "dashed") +
		labs(x = "RNA1", y = "RNA2", color = "Gamma") +
		scale_color_discrete(labels = c("0.00",
																		"0.01", "0.05", "0.10", "0.15",
																		"0.20", "0.25", "0.30", "0.35",
																		"0.40", "0.45", "0.50", "0.53")) +
		scale_x_continuous(limits = c(NA, limits[1]), breaks = seq(0, 1000, 100)) +
		scale_y_continuous(limits = c(NA, limits[2]), breaks = seq(0, 1000, 100)) +
		theme_bw()

	if (zoom.box) {
		gg <- gg +
			annotate("rect", xmin = 0, xmax = limits[1], ymin = 0, ymax = limits[2],
							colour = "springgreen3", alpha = 0.2) +
			scale_x_continuous(limits = c(NA, 1000), breaks = seq(0, 1000, 100)) +
			scale_y_continuous(limits = c(NA, 1000), breaks = seq(0, 1000, 100))
	}
	return(gg)
}

# plot_6eq_separatrix(df.gg.init, df.gg.sep, drift = c(5, -3))
plot_6eq_separatrix(df.gg.init, df.gg.sep, drift = c(5, -3), limits = c(500,300), zoom.box = T) #+ theme_bw(base_family = "LM Roman 10") + alfR::ggexport("gg_6eq_separatrix_point_full.pdf", size = c(6, 4), font = "lmodern")
plot_6eq_separatrix(df.gg.init, df.gg.sep, drift = c(3, -1), limits = c(500,300)) #+ theme_bw(base_family = "LM Roman 10") + alfR::ggexport("gg_6eq_separatrix_point_zoom.pdf", size = c(6, 4), font = "lmodern")


gg_6eq_separatrix_point + theme_bw(base_family = "LM Roman 10") + alfR::ggexport("gg_6eq_separatrix_point_full.pdf", size = c(6, 4), font = "lmodern")
# saveRDS(gg_6eq_separatrix_point, paste0("objects/", "gg_2eq_separatrix_point_parR2_10.rds"))
