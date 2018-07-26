library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions


# Parameters
kappa <- 1
alpha <- 1
sigma <- 0.1

# Number of interations
n.iter <- 2000
gamma.parts <- 20
R.parts <- 1000
p.parts <- 1000

params <- c(n.iter,
						0, 0,
						kappa, alpha, 0, sigma,
						gamma.parts, R.parts, p.parts)

system("gcc -Ofast -lm 2eq_rk4_bipartite_var_gamma.c -o compiled/rk4-gamma-var")
alfR::lok_regar(system(paste("./compiled/rk4-gamma-var", paste(params, collapse = ' '), "> data/results-gamma-var.csv")))

data2eq.gm20 <- data.table::fread("data/results-gamma-var.csv") %>%
	dplyr::rename(gamma = V1, R.init = V2, p.init = V3) %>%
	dplyr::arrange(desc(gamma)) %>%
	dplyr::filter(
		gamma <= 0.55,
		gamma != 0
	)




gg.2eq.gamma20 <- ggplot(data2eq.gm20) +
	geom_tile(aes(x = R.init, y = p.init, fill = (gamma))) +
	scale_fill_gradient(low = "seagreen4", high = "darkseagreen1",
											breaks = seq(0.1, 0.5, 0.05),
											guide = "legend") +
	labs(fill = "Gamma") +
	theme_bw()

gg.2eq.gamma20

# plot.2eq.gm10 <- readRDS("objects/gg_2eq_gamma10.rds")

saveRDS(data2eq.gm20, "objects/data2eq.gm20.rds")

saveRDS(gg.2eq.gamma20, "objects/gg.2eq.gamma20.rds")
