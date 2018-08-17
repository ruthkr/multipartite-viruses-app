library(tidyverse)

source("2eq_find_first_point_to_P2.R")


# Functions ---------------------------------------------------------------

get_det_data <- function(R1.init, p.init, gamma, scale = 1000) {
	# R.init <- 0.408
	# p.init <- 0

	# Parameters
	kappa <- 1
	omega <- 1
	# gamma <- 0.5
	sigma <- 0.1

	# Number of interations
	n.iter <- 1500

	# Flags
	# fixed.point.text <- "plus"
	# fixed.point.select <- paste0("second.", fixed.point.text)

	params <- c(n.iter,
							R1.init, p.init,
							kappa, omega, gamma, sigma)

	#Compile
	# system("gcc -Ofast -lm 2eq_rk4_bipartite_test.c -o compiled/rk4-test")
	system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

	# Read data
	data.2eq <- data.table::fread("data/results-test.csv")

	return(data.2eq * scale)
}

get_sto_data <- function(R1.init, p.init, gamma, time.max = 100) {
	Sys.sleep(1)
	# Initial conditions
	# R1.init <- 408
	# p.init <- 0

	# Parameters
	kappa <- 1
	omega <- 1
	# gamma <- 0.5
	sigma <- 0.1

	# Number maximum time
	# time.max <- 100

	params <- c(time.max,
							R1.init, p.init,
							kappa, omega, gamma, sigma)

	#Compile
	# system("gcc -Ofast -lm 2eq_ssa_new.c -o compiled/2eq-ssa")
	system(paste("./compiled/2eq-ssa", paste(params, collapse = ' '), "> data/results_ssa_c.csv"))

	# Read data
	data.2eq.ssa <- data.table::fread("data/results_ssa_c.csv")

	return(data.2eq.ssa)
}


# Plot --------------------------------------------------------------------

scale <- 1000

nullclines.R <- function(kappa, gamma, num) {
	p = 0;
	R = numeric(num+1);
	p.step = 1/num;
	for(i in 0:num+1){
		R[i] = (kappa * p - gamma) / (kappa * p)
		p = p + p.step;
	}
	return(R)
}

nullclines.p <- function(omega, sigma, num) {
	R = 0;
	p = numeric(num+1);
	R.step = 1/num;
	for(i in 0:num+1) {
		p[i] <- (omega * R) / (omega * R + sigma)
		R = R + R.step;
	}
	return(p)
}

plot_trajectory <- function(kappa, omega, gamma, sigma, step = 1000, scale = 1000, show.sto = F) {
	time <- seq(0, 1, 1/step) * scale

	# Nullclines
	R <- nullclines.R(kappa, gamma, step) * scale
	p <- nullclines.p(omega, sigma, step) * scale
	df.nullclines <- data.frame(time, R, p)

	a <- kappa * (sigma/omega + 1)
	b <- -(kappa + gamma)
	c <- gamma

	p.plus <- ((-b + sqrt(b^2 - 4*a*c)) / (2*a))
	R.plus <- (1 - gamma/(kappa * p.plus))
	p.minus <- ((-b - sqrt(b^2 - 4*a*c)) / (2*a))
	R.minus <- (1 - gamma/(kappa * p.minus))

	p.plus <- p.plus * scale
	R.plus <- R.plus * scale
	p.minus <- p.minus * scale
	R.minus <- R.minus * scale

	# Make plot
	gg <- ggplot(df.nullclines) +
		# Nullclines
		geom_line(data = df.nullclines %>% filter(R >= 0), aes(x = R , y = time, colour = "RNA", alpha = "RNA")) +
		geom_line(aes(x = time , y = p, colour = "viral.replicase", alpha = "viral.replicase")) +
		# Separatrix
		geom_line(data = df.clean %>%
								filter(gamma == 0.5) %>%
								mutate(R1.init = R1.init * scale,
											 p.init = p.init * scale),
							aes(x = (R1.init - 1), y = p.init, colour= "separatrix", size = "separatrix"), linetype = "longdash")
		# Paths for different initial conditions
		if (show.sto) {
			# Stochastic
			gg <- gg +
				geom_path(data = get_sto_data(100, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.100", linetype = "sto")) +
				geom_path(data = get_sto_data(250, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.250", linetype = "sto")) +
				geom_path(data = get_sto_data(408, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.408", linetype = "sto")) +
				geom_path(data = get_sto_data(550, 0, 0.5,  25), aes(x = V2 + 1, y = V3, colour = "path.550", linetype = "sto")) +
				geom_path(data = get_sto_data(700,  0, 0.5, 15), aes(x = V2 + 1, y = V3, colour = "path.700", linetype = "sto")) +
				geom_path(data = get_sto_data(850,  0, 0.5, 15), aes(x = V2 + 1, y = V3, colour = "path.850", linetype = "sto")) +
				geom_path(data = get_sto_data(1000, 0, 0.5, 15), aes(x = V2 + 1, y = V3, colour = "path.1000", linetype = "sto"))
		}
		# Deterministic
	gg <- gg +
		geom_path(data = get_det_data(0.100, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.100", linetype = "det")) +
		geom_path(data = get_det_data(0.250, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.250", linetype = "det")) +
		geom_path(data = get_det_data(0.408, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.408", linetype = "det")) +
		geom_path(data = get_det_data(0.550, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.550", linetype = "det")) +
		geom_path(data = get_det_data(0.700, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.700", linetype = "det")) +
		geom_path(data = get_det_data(0.850, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.850", linetype = "det")) +
		geom_path(data = get_det_data(1.000, 0, 0.5), aes(x = V2 + 1, y = V3, colour = "path.1000", linetype = "det")) +
		# Initial conditions
		geom_point(aes(x = 100, y = 0, colour = "path.100")) +
		geom_point(aes(x = 250, y = 0, colour = "path.250")) +
		geom_point(aes(x = 408, y = 0, colour = "path.408")) +
		geom_point(aes(x = 550, y = 0, colour = "path.550")) +
		geom_point(aes(x = 700, y = 0, colour = "path.700")) +
		geom_point(aes(x = 850, y = 0, colour = "path.850")) +
		geom_point(aes(x = 1000, y = 0, colour = "path.1000")) +
		# Fixed points
		geom_point(aes(x = R.minus, y = p.minus), colour = "black") +
		geom_point(aes(x = R.plus, y = p.plus), colour = "black") +
		# Extra ggplot stuff
		labs(x = "RNA1", y = "Replicase",
				 alpha = "Nullclines", color = "Initial conditions",
				 linetype = "Simulation type", size = "") +
		scale_size_manual(values = c("separatrix" = 0.5), labels = c("separatrix" = "Separatrix")) +
		scale_alpha_manual(values = c("RNA" = 1,
																	"viral.replicase" = 1
																	),
											labels = c("RNA" = "RNA1",
																 "viral.replicase" = "Replicase"
											)) +
		scale_color_manual(values = c("RNA" = "red3",
																	"viral.replicase" = "springgreen4",
																	"separatrix" = "black",
																	# Paths
																	"path.100" = "palevioletred1",
																	"path.250" = "orchid2",
																	"path.408" = "purple3",
																	"path.550" = "slateblue4",
																	"path.700" = "dodgerblue3",
																	"path.850" = "royalblue1",
																	"path.1000" = "steelblue1"
																	),
											breaks = c("path.100",
																 "path.250",
																 "path.408",
																 "path.550",
																 "path.700",
																 "path.850",
																 "path.1000"
											),
											labels = c("path.100"  = expression(paste(R["1"], " = 100  ")),
																 "path.250"  = expression(paste(R["1"], " = 250  ")),
																 "path.408"  = expression(paste(R["1"], " = 408  ")),
																 "path.550"  = expression(paste(R["1"], " = 550  ")),
																 "path.700"  = expression(paste(R["1"], " = 700  ")),
																 "path.850"  = expression(paste(R["1"], " = 850  ")),
																 "path.1000" = expression(paste(R["1"], " = 1000"))
											)) +
		scale_linetype_manual(values = c("sto" = "solid",
																		 "det" = "dashed"
																		 ),
													labels = c("sto" = "Stochastic",
																		 "det" = "Deterministic"
																		 )) +
		scale_x_continuous(breaks = seq(0, 1000, 100)) +
		guides(
			size = guide_legend(order = 1),
			alpha = guide_legend(order = 2,
													 override.aes = list(colour = c("RNA" = "red3",
													 															 "viral.replicase" = "springgreen4"))),
			color = guide_legend(order = 3,
													 override.aes = list(shape = c(NA)))
			) +
		theme_bw()

	return(gg)
}

# (gg.trajs.det <- plot_trajectory(1, 1, 0.5, 0.1))
# (gg.trajs.det.sto <- plot_trajectory(1, 1, 0.5, 0.1, show.sto = TRUE))

# saveRDS(gg.trajs.det, "objects/2eq_nullcline_det.rds")
# saveRDS(gg.trajs.det.sto, "objects/2eq_nullcline_det_sto_final.rds")


# Read plots --------------------------------------------------------------

gg.trajs.det <- readRDS("objects/2eq_nullcline_det.rds")
gg.trajs.det.sto <- readRDS("objects/2eq_nullcline_det_sto_final.rds")

gg.trajs.det #+ theme(text = element_text(family = "LM Roman 10")) + ggsave(filename = "2eq_nullcline_det.pdf", width = 6, height = 4.5, dpi = 96, device = "pdf")
gg.trajs.det.sto #+ theme(text = element_text(family = "LM Roman 10")) + ggsave(filename = "2eq_nullcline_det_sto.pdf", width = 6, height = 4.5, dpi = 96, device = "pdf")
