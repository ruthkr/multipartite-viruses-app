library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions


# Parameters
gamma <- 0.1
alpha <- 1
sigma <- 0.1

# Number of interations
n.iter <- 2000
kappa.parts <- 20
R.parts <- 1000
p.parts <- 1000

params <- c(n.iter,
						0, 0,
						0, alpha, gamma, sigma,
						kappa.parts, R.parts, p.parts)

system("gcc -Ofast -lm 2eq_rk4_bipartite_var_kappa.c -o compiled/rk4-kappa-var")
alfR::lok_regar(system(paste("./compiled/rk4-kappa-var", paste(params, collapse = ' '), "> data/results-kappa-var.csv")))

 %>%  <- data.table::fread("data/results-kappa-var.csv") %>%
	dplyr::rename(kappa = V1, R.init = V2, p.init = V3) %>%
	dplyr::arrange(desc(kappa)) %>%
	dplyr::filter(
		# kappa <= 1,
		kappa > 0.80
	)


# data2eq.kp20 <- readRDS("objects/data2eq.kp20.rds")

gg.2eq.kappa20 <- ggplot(data2eq.kp20) +
	geom_tile(aes(x = R.init, y = p.init, fill = (kappa))) +
	scale_fill_gradient(low = "darkseagreen1", high = "seagreen4",
											breaks = seq(0.9, 1, 0.1),
											guide = "legend") +
	labs(fill = "Kappa") +
	theme_bw()

gg.2eq.kappa20

# plot.2eq.gm10 <- readRDS("objects/gg_2eq_gamma10.rds")

# saveRDS(data2eq.gm20, "objects/data2eq.gm20.rds")

saveRDS(gg.2eq.gamma20, "objects/plot_2eq_var_gamma20.rds")
