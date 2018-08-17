library(tidyverse)

# Plot trajectories from the results of deterministic simualtion

plot_traj_6eq <- function(data.sto, data.det, extra.data = NULL){
	# Base plot
	gg <- ggplot(data.det) +
		aes(x = V2, y = V3) +
		geom_path(aes(colour = "deterministic")) +
		geom_point(data = data.sto[1,], aes(shape = "initial"), size = 3.5) +
		geom_point(data = data.det[nrow(data.det),], aes(colour = "deterministic", shape = "fixed"), size = 2.5) +
		labs(x = "RNA1", y = "RNA2", color = "Simulation Path", shape = "Point Type") +
		scale_color_manual(values = c("others" = "springgreen4",
																	"stochastics" = "mediumorchid",
																	"deterministic" = "steelblue1")) +
		scale_shape_manual(values = c("fixed" = 16,
																	"initial" = 18)) +
		guides(color = guide_legend(order = 1, override.aes = list(shape = c(NA,NA)))) +
		theme_bw()

	# Extra plots
	if (length(extra.data) != 0) {
		for (i in 1:length(extra.data)) {
			gg <- gg +
				geom_line(data = extra.data[[i]], aes(colour = "deterministic"), linetype = 2) +
				geom_point(data = extra.data[[i]][1,], aes(shape = "initial"), size = 3.5) +
				geom_point(data = extra.data[[i]][nrow(extra.data[[i]]),], aes(colour = "deterministic", shape = "fixed"), size = 2.5)
		}
	}

	# Add stochastic path
	gg <- gg +
		geom_path(data = data.sto, aes(colour = "stochastics")) +
	geom_point(data = data.sto[nrow(data.sto),], aes(colour = "stochastics", shape = "fixed"), size = 2.5)

	return(gg)
}

# Example of usage

# For R1=100, R2=104
# data.6eq.ssa.100.104 <- readRDS("objects/data.6eq.ssa.100.104.rds")
# data.6eq.det.01.0104 <- readRDS("objects/data.det.6eq.01.0104.rds")
# plot_traj_6eq(data.6eq.ssa.500.200, data.6eq.det.05.02)
# plot_traj_6eq(data.6eq.ssa.100.200, data.6eq.det.01.02)
# plot_traj_6eq(
# 	data.6eq.ssa.100.104, data.6eq.det.01.0104,
# 	list(data.det.6eq.005.02, data.det.6eq.01.02, data.det.6eq.015.02)
# )

# For R1=400, R2=500
# vars.factor <- 1000
# data.6eq.ssa.400.500 <- data.table::fread("data/results_6eq_ssa_c_400_500_gamma02.csv")
# data.6eq.det.400.500 <- data.table::fread("data/results-test_040_050_gamma02.csv")
# data.6eq.det.400.500 <- data.6eq.det.400.500 * vars.factor
# plot.trajectories.2D.400.500 <- plot_traj_6eq(
# 	data.6eq.ssa.400.500, data.6eq.det.400.500,
# 	list(data.det.6eq.039.05, data.det.6eq.025.05, data.det.6eq.01.05)
# )

# For R1=400, R2=500
# vars.factor <- 1000
# data.6eq.ssa.400.250 <- data.table::fread("data/results_6eq_ssa_c_400_250_gamma03.csv")
# data.6eq.det.400.250 <- data.table::fread("data/results-test_040_025_gamma03.csv")
# data.6eq.det.400.250 <- data.6eq.det.400.250 * vars.factor
# plot.trajectories.2D.400.250 <- plot_traj_6eq(
# 	data.6eq.ssa.400.250, data.6eq.det.400.250,
# 	list(data.det.6eq.038.025, data.det.6eq.01.025, data.det.6eq.03.025)
# )
#
# plot.trajectories.2D.400.250
