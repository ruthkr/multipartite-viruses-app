library(tidyverse)

# Number maximum time
time.max <- 2000

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 1000
R2.init <- 800
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

gg.data.6eq.ssa  <- ggplot(data.6eq.ssa) +
	# Stochastic
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "R2")) +
	geom_line(aes(x = V1, y = V4, colour = "p")) +
	geom_line(aes(x = V1, y = V5, colour = "S")) +
	geom_line(aes(x = V1, y = V6, colour = "v1")) +
	geom_line(aes(x = V1, y = V7, colour = "v2")) +
	# # Deterministic
	# geom_line(data = data.6eq.det, aes(x = V1, y = V2, colour = "R1"), linetype = "dashed") +
	# geom_line(data = data.6eq.det, aes(x = V1, y = V3, colour = "R2"), linetype = "dashed") +
	# geom_line(data = data.6eq.det, aes(x = V1, y = V4, colour = "p"), linetype = "dashed") +
	# geom_line(data = data.6eq.det, aes(x = V1, y = V5, colour = "S"), linetype = "dashed") +
	# geom_line(data = data.6eq.det, aes(x = V1, y = V6, colour = "v1"), linetype = "dashed") +
	# geom_line(data = data.6eq.det, aes(x = V1, y = V7, colour = "v2"), linetype = "dashed") +
	labs(x = "Time (Stochastic)", y = "Variables", color = "Variable") +
	# scale_x_continuous(sec.axis = sec_axis(~.*10, name = "Time (Deterministic)"), limits = c(0, 90)) +
	scale_color_manual(values = c("R1" = "red3",        "R2" = "royalblue",
																"p" = "springgreen4", "S" = "purple3",
																"v1" = "black", "v2" = "orange")) +
	theme_bw()

gg.data.6eq.ssa

tail(data.6eq.ssa)

# #  Save Data
# saveRDS(data.6eq.det, "objects/data.6eq.det.rds")
# saveRDS(data.6eq.ssa, "objects/data.6eq.ssa.500.200.exit.rds")
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
