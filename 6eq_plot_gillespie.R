library(tidyverse)

# Number maximum time
time.max <- 100

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 160
R2.init <- 200
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

params <- c(time.max,
						R1.init, R2.init, p.init, s.init, v1.init, v2.init,
						kappa, omega, gamma, sigma, epsilon, delta,
						current.capacity)

#Compile
system("gcc -Ofast -lm 6eq_ssa.c -o compiled/6eq-ssa")
alfR::lok_regar(system(paste("./compiled/6eq-ssa", paste(params, collapse = ' '), "> data/results_6eq_ssa_c.csv")))


# Read data
data.6eq.ssa <- data.table::fread("data/results_6eq_ssa_c.csv")
saveRDS(data.6eq.ssa, "data/6eq.sto.R1160.R2200.gamma02.rds")

gg.data.6eq.ssa  <- ggplot(data.6eq.ssa) +
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "R2")) +
	geom_line(aes(x = V1, y = V4, colour = "p")) +
	geom_line(aes(x = V1, y = V5, colour = "s")) +
	geom_line(aes(x = V1, y = V6, colour = "v1")) +
	geom_line(aes(x = V1, y = V7, colour = "v2")) +
	labs(x = "Time", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R1" = "red3", "R2" = "royalblue",
																"p" = "springgreen4", "s" = "purple3",
																"v1" = "black", "v2" = "orange"),
										 breaks = c("R1", "R2", "p", "s", "v1", "v2"),
										 labels = c("R1", "R2", "p", "s", "v1", "v2")) +
	theme_bw()

gg.data.6eq.ssa

# tail(data.6eq.ssa)

# #  Save Data


saveRDS(gg.data.6eq.ssa, "objects/gg.6eq.sto.R1160.R2200.gamma02.rds")

# saveRDS(data.6eq.det, "objects/data.6eq.det.rds")
# saveRDS(gg.data.6eq.ssa, "objects/gg.data.6eq.ssa.500.200.exit.rds")
# gg.data.6eq.ssa <- readRDS("objects/gg.data.6eq.ssa.rds")



# gg.data.2eq.ssa.420 + ggsave(filename = "ssa_420.png", width = 6, height = 3, dpi = 300, device = "png")

# Trajectories

# ggplot(data.2eq.ssa.401.3) +
# 	geom_path(aes(x = V2, y = V3, colour = "R")) +
# 	labs(x = "RNA1", y = "p", color = "Variable") +
# 	scale_color_manual(values = c("R" = "red3",
# 																"p" = "springgreen4")) +
# 	theme_bw()
