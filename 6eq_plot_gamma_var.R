library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

p.init <- 0
s.init <- 0
v1.init <- 0
v2.init <- 0

# Parameters
kappa1 <- 1
kappa2 <- 1

alpha <- 1
beta <- 1

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 2500
gamma.parts <- 10
R.parts <- 500

params <- c(n.iter,
						0, 0, p.init, s.init,
						v1.init, v2.init,
						kappa1, kappa2, alpha, beta, 0, 0, sigma1, sigma2,
						epsilon1, epsilon2, delta1, delta2,
						gamma.parts, R.parts)

system("gcc -Ofast -lm 6eq_rk4_bipartite_var_gamma.c -o compiled/rk4-gamma-var")
alfR::lok_regar(system(paste("./compiled/rk4-gamma-var", paste(params, collapse = ' '), "> data/results-gamma-var.csv")))

data_6eq_var_gamma10_until06 <- data.table::fread("data/results-gamma-var.csv") %>%
	dplyr::rename(gamma = V1, R1.init = V2, R2.init = V3) %>%
	dplyr::arrange(desc(gamma)) %>%
	dplyr::filter(
		gamma <= 0.6,
		gamma != 0
		)

plot_6eq_var_gamma10 <- ggplot(data_6eq_var_gamma10) +
	geom_tile(aes(x = R1.init, y = R2.init, fill = (gamma))) +
	scale_fill_gradient(low = "seagreen4", high = "darkseagreen1",
											breaks = seq(0.1, 0.5, 0.1),
											guide = "legend") +
	labs(fill = "Gamma", x = "Initial Value of RNA1", y = "Initial Value of RNA2") +
	theme_bw()

saveRDS(plot_6eq_var_gamma10_until06, "objects/plot_6eq_var_gamma10_until06.rds")


