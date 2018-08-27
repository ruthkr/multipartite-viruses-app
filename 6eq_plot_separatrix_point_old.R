library(tidyverse)

source("6eq_find_first_point_to_P3.R")

scale <- 1000

get_det_data <- function(R1.init, R2.init, gamma, scale = 1000) {

	p.init <- 0
	s.init <- 0
	v1.init <- 0
	v2.init <- 0

	# Parameters
	kappa1 <- 1
	kappa2 <- 1

	alpha <- 1
	beta <- 1
	omega <- 1

	sigma1 <- 0.1
	sigma2 <- 0.1

	epsilon1 <- 0.1
	epsilon2 <- 0.1

	delta1 <- 0.1
	delta2 <- 0.1

	# Number of interations
	options(scipen = 999)
	n.iter <- 1500

	# Flags
	# fixed.point.select <- "third.plus"

	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							v1.init, v2.init,
							kappa1, kappa2, alpha, beta, gamma, gamma, sigma1, sigma2,
							epsilon1, epsilon2, delta1, delta2)

	#Compile
	system("gcc -Ofast -lm 6eq_rk4_bipartite_test.c -o compiled/rk4-test")
	system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

	# Read data
	data.6eq <- data.table::fread("data/results-test.csv")

	return(data.6eq * scale)
}

df.gg.sep <- df.clean.sep %>%
	filter(gamma == 0.2) %>%
	mutate(R1.init = R1.init * scale,
				 R2.init = R2.init * scale) %>%
	group_by(gamma, R1.init) %>%
	summarise(R2.init = min(R2.init))


df.gg.init <- df.clean.init %>%
	filter(gamma == 0.2) %>% filter(R2.init %in% c(0.2, 0.4, 0.6, 0.8, 1))%>%
	mutate(R1.init = R1.init * scale,
				 R2.init = R2.init * scale)

data.6eq.sto.R1160.R2200.gamma02 <-  readRDS("data/6eq.sto.R1160.R2200.gamma02.rds")

gg_6eq_separatrix_point <- ggplot(df.gg.init) +
	geom_line(data = df.gg.sep, aes(x = (R1.init - 1), y = R2.init, colour= "separatrix"), linetype = "longdash") +
	labs(x = "RNA1", y = "RNA2", color = "Simulation Type") +
	scale_x_continuous(breaks = seq(0,1000,100)) +
	scale_y_continuous(breaks = seq(0,1000,100)) +
	scale_color_manual(values = c("stochastic" = "mediumorchid",
																"deterministic" = "steelblue1",
																"separatrix" = "black")) +
	scale_linetype_manual(values = c("sto" = "solid",
																	"det" = "dashed"
																	),
												guide = FALSE) +
	guides(color = guide_legend(override.aes = list(shape = c(NA,NA,NA)))) +
	theme_bw()

gg <- gg_6eq_separatrix_point +
	geom_path(data = get_det_data(0.2, 1, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(0.4, 1, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(0.6, 1, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(1, 0.8, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(1, 0.6, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(1, 0.4, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(1, 0.2, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	geom_path(data = get_det_data(1, 1, 0.2), aes(x = V2 + 1, y = V3, colour = "deterministic", linetype = "det")) +
	# Initial conditions
	geom_point(aes(x = 200, y = 1000, colour = "deterministic")) +
	geom_point(aes(x = 400, y = 1000, colour = "deterministic")) +
	geom_point(aes(x = 600, y = 1000, colour = "deterministic")) +
	geom_point(aes(x = 1000, y = 1000, colour = "deterministic")) +
	geom_point(aes(x = 1000, y = 800, colour = "deterministic")) +
	geom_point(aes(x = 1000, y = 600, colour = "deterministic")) +
	geom_point(aes(x = 1000, y = 400, colour = "deterministic")) +
	geom_point(aes(x = 1000, y = 200, colour = "deterministic")) +
	geom_point(aes(x = 160, y = 200, colour = "deterministic")) +
	# Add stochastic path
	geom_path(data = data.6eq.sto.R1160.R2200.gamma02, aes(x = V2 + 2, y = V3, colour = "stochastic", linetype = "sto")) +
	#Near separatrix
	geom_path(data = get_det_data(0.16, 0.2, 0.2), aes(x = V2 + 2, y = V3, colour = "deterministic", linetype = "det"))


gg #+ alfR::ggexport("gg_6eq_separatrix_point_R1160_R2200.pdf", font = "lmodern")
# saveRDS(gg_6eq_separatrix_point, paste0("objects/", "gg_2eq_separatrix_point_parR2_10.rds"))

# %in% c(0.00, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53
