# This Script is used to plot the Deterministic and Stochastic Simulation in the same time

library(tidyverse)

# Deterministic Simulation
# Parameters --------------------------------------------------------------
# Initial conditions

R1.init.det <- 0.1
R2.init.det <- 0.104
p.init.det <- 0
s.init.det <- 0
v1.init.det <- 0
v2.init.det <- 0

# Parameters
kappa1.det <- 1
kappa2.det <- 1

alpha.det <- 1
beta.det <- 1

gamma1.det <- 0.20
gamma2.det <- 0.20

sigma1.det <- 0.1
sigma2.det <- 0.1

epsilon1.det <- 0.1
epsilon2.det <- 0.1

delta1.det <- 0.1
delta2.det <- 0.1

# Number of interations
n.iter <- 7500

params.det <- c(n.iter,
						R1.init.det, R2.init.det, p.init.det, s.init.det,
						v1.init.det, v2.init.det,
						kappa1.det, kappa2.det, alpha.det, beta.det, gamma1.det, gamma2.det, sigma1.det, sigma2.det,
						epsilon1.det, epsilon2.det, delta1.det, delta2.det)

#Compile
system("gcc -Ofast -lm 6eq_rk4_bipartite_test.c -o compiled/rk4-test")
system(paste("./compiled/rk4-test", paste(params.det, collapse = ' '), "> data/results-test.csv"))

# Read data
data.det <- data.table::fread("data/results-test.csv")
data.det.6eq <- data.det * 1000

# Rename data file
det <- assign(paste("data.det.6eq", gsub("[.]", "", R1.init.det), gsub("[.]", "", R2.init.det), sep = "."), data.det.6eq)

# Stochastics Simulation
# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 100
R2.init <- 104
p.init <- 0
s.init <- 0
v1.init <- 0
v2.init <- 0

# Parameters
kappa <- 1
omega <- 1
gamma <- 0.2
sigma <- 0.1
epsilon <- 0.1
delta <- 0.1

# Current Capacity
current.capacity <- 1000

params.sto <- c(time.max,
						R1.init, R2.init, p.init, s.init, v1.init, v2.init,
						kappa, omega, gamma, sigma, epsilon, delta,
						current.capacity)

#Compile
system("gcc -Ofast -lm 6eq_ssa.c -o compiled/6eq-ssa")
system(paste("./compiled/6eq-ssa", paste(params.sto, collapse = ' '), "> data/results_6eq_ssa_c.csv"))

# Read data
data.ssa.6eq <- data.table::fread("data/results_6eq_ssa_c.csv")

# Rename data file name
ssa <- assign(paste("data.ssa.6eq", R1.init, R2.init, sep = "."), data.ssa.6eq)


# Plot both Stochastics and Deterministic Result Simulation
gg.6eq.ssa.det  <- ggplot(ssa) +
	# Stochastic
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "R2")) +
	geom_line(aes(x = V1, y = V4, colour = "p")) +
	geom_line(aes(x = V1, y = V5, colour = "S")) +
	geom_line(aes(x = V1, y = V6, colour = "v1")) +
	geom_line(aes(x = V1, y = V7, colour = "v2")) +
	# Deterministic
	geom_line(data = det, aes(x = V1, y = V2, colour = "R1"), linetype = "dashed") +
	geom_line(data = det, aes(x = V1, y = V3, colour = "R2"), linetype = "dashed") +
	geom_line(data = det, aes(x = V1, y = V4, colour = "p"), linetype = "dashed") +
	geom_line(data = det, aes(x = V1, y = V5, colour = "S"), linetype = "dashed") +
	geom_line(data = det, aes(x = V1, y = V6, colour = "v1"), linetype = "dashed") +
	geom_line(data = det, aes(x = V1, y = V7, colour = "v2"), linetype = "dashed") +
	labs(x = "Time (Stochastic)", y = "Variables", color = "Variable") +
	scale_x_continuous(sec.axis = sec_axis(~.*10, name = "Time (Deterministic)"), limits = c(0, 100)) +
	scale_color_manual(values = c("R1" = "red3",        "R2" = "royalblue",
																"p" = "springgreen4", "S" = "purple3",
																"v1" = "black", "v2" = "orange")) +
	theme_bw()

assign(paste("gg.6eq.ssa.det", R1.init, R2.init, sep = "."), gg.6eq.ssa.det)
