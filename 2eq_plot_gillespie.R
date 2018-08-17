library(tidyverse)
#
#
# # Simulation Using Python
# data.ssa <- data.table::fread("data/results-ssa.csv")
#
# ggplot(data.ssa) +
# 	geom_line(aes(x = V1, y = V2, colour = "R")) +
# 	geom_line(aes(x = V1, y = V3, colour = "p")) +
# 	labs(x = "Time", y = "Variables", color = "Variable") +
# 	scale_color_manual(values = c("R" = "red3",
# 																"p" = "springgreen4")) +
# 	theme_bw()

# Simulation Using C
# Parameters --------------------------------------------------------------
# Initial conditions
R.init <- 550
p.init <- 0

# Parameters
kappa <- 1
omega <- 1
gamma <- 0.5
sigma <- 0.1

# Number maximum time
time.max <- 1000

params <- c(time.max,
						R.init, p.init,
						kappa, omega, gamma, sigma)

#Compile
system("gcc -Ofast -lm 2eq_ssa_new.c -o compiled/2eq-ssa")
system(paste("./compiled/2eq-ssa", paste(params, collapse = ' '), "> data/results_ssa_c.csv"))

# Read data
data.2eq.ssa <- data.table::fread("data/results_ssa_c.csv")

# gg.data.2eq.ssa.401.2 <- ggplot(
# 	rbind(
# 		data.2eq.ssa.401.2[1:200],
# 		sample_n(data.2eq.ssa.401.2[201:(nrow(data.2eq.ssa.401.2) - 200)], 8000),
# 		data.2eq.ssa.401.2[(nrow(data.2eq.ssa.401.2) - 199):nrow(data.2eq.ssa.401.2)]
# 		)
# 	) +
# 	geom_line(aes(x = V1, y = V2, colour = "R")) +
# 	geom_line(aes(x = V1, y = V3, colour = "p")) +
# 	labs(x = "Time", y = "Variables", color = "Variable") +
# 	scale_color_manual(values = c("R" = "red3",
# 																"p" = "springgreen4")) +
# 	theme_bw()

gg.ssa <- ggplot(data.2eq.ssa) +
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "p")) +
	labs(x = "Time", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R1" = "red3",
																"p" = "springgreen4")) +
	theme_bw()

gg.ssa

saveRDS(gg.ssa,  paste0("objects/","gg.2eq.ssa.R550.p0.gamma05.rds"))

# saveRDS(gg.data.2eq.ssa.R400.p10.gamma05, paste0("objects/", "gg.data.2eq.ssa.R400.p10.gamma05.rds"))



# saveRDS(gg.data.2eq.ssa.420, "objects/gg.data.2eq.ssa.420.rds")
# gg.data.2eq.ssa.420 <- readRDS("objects/gg.data.2eq.ssa.420.rds")
#
# gg.data.2eq.ssa.420 + ggsave(filename = "ssa_420.png", width = 6, height = 3, dpi = 300, device = "png")

# Trajectories

# ggplot(data.2eq.ssa.401.3) +
# 	geom_path(aes(x = V2, y = V3, colour = "R")) +
# 	labs(x = "RNA1", y = "p", color = "Variable") +
# 	scale_color_manual(values = c("R" = "red3",
# 																"p" = "springgreen4")) +
# 	theme_bw()
