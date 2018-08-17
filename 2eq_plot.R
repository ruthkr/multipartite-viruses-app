library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 0.408
p.init <- 0

# Parameters
kappa <- 1
omega <- 1
gamma <- 0.5
sigma <- 0.1

# Number of interations
n.iter <- 1500

# Flags
fixed.point.text <- "plus"
fixed.point.select <- paste0("second.", fixed.point.text)

params <- c(n.iter,
						R1.init, p.init,
						kappa, omega, gamma, sigma)

#Compile
system("gcc -Ofast -lm 2eq_rk4_bipartite_test.c -o compiled/rk4-test")
system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

# Read data
data.2eq <- data.table::fread("data/results-test.csv")

# Scale
scale <- 1000

# # Create plot
# gg <- ggplot(data.2eq) +
# 	geom_line(aes(x = V1, y = V2, colour = "RNA")) +
# 	geom_line(aes(x = V1, y = V3, colour = "viral.replicase")) +
# 	labs(x = "Time", y = "Variables", color = "Variable") +
# 	scale_color_manual(values = c("RNA" = "red3",
# 																"viral.replicase" = "springgreen4")) +
# 	theme_bw()
#
# gg


# Second Fixed Point
a <- kappa * (sigma/omega + 1)
b <- -(kappa + gamma)
c <- gamma

p.plus <- (-b + sqrt(b^2 - 4*a*c)) / (2*a)
R.plus <- 1 - gamma/(kappa * p.plus)

p.minus <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)
R.minus <- 1 - gamma/(kappa * p.minus)


# if (fixed.point.select == "second.plus") {
#   gg.second.plus <- gg +
# 		geom_hline( yintercept = p.plus, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R.plus, color = "red3", linetype="dotted" )
# } else if (fixed.point.select == "second.minus") {
# 	gg.second.minus <- gg +
# 		geom_hline( yintercept = p.minus, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R.minus, color = "red3", linetype="dotted" )
# }
#
# gg.second.plus

# get(paste0("gg.second.", fixed.point.text))
#
# saveRDS(get(paste0("gg.second.", fixed.point.text)),
# 				paste0("objects/", "timeseries_2eq_",
# 							 "R1", gsub("[.]", "", R1.init), "_",
# 							 "p", gsub("[.]", "", p.init), "_",
# 							 "gamma", gsub("[.]", "", gamma), "_",
# 							 "p2", fixed.point.text,
# 							 ".rds"))

# ggplot(data.2eq) +
# 	geom_path(aes(x = V2, y = V3, colour = "R")) +
# 	labs(x = "Time", y = "Variables", color = "Variable") +
# 	scale_color_manual(values = c("R" = "red3",
# 																"p" = "springgreen4")) +
# 	theme_bw()

# Create plot with the scale
p.plus.scaled <- p.plus * scale
R.plus.scaled <- R.plus * scale
p.minus.scaled <- p.minus * scale
R.minus.scaled <- R.minus * scale

gg <- ggplot(data.2eq %>% mutate(V2 = V2 * scale, V3 = V3 * scale)) +
	geom_line(aes(x = V1, y = V2, colour = "RNA")) +
	geom_line(aes(x = V1, y = V3, colour = "viral.replicase")) +
	labs(x = "Time", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("RNA" = "red3",
																"viral.replicase" = "springgreen4")) +
	theme_bw()

if (fixed.point.select == "second.plus") {
  gg.second.plus <- gg +
		geom_hline( yintercept = p.plus.scaled, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R.plus.scaled, color = "red3", linetype="dotted" )
} else if (fixed.point.select == "second.minus") {
	gg.second.minus <- gg +
		geom_hline( yintercept = p.minus.scaled, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R.minus.scaled, color = "red3", linetype="dotted" )
}

gg.second.plus


saveRDS(gg.second.plus, paste0("objects/", "gg.2eq.det.R480.p0.gamma05.rds"))


# gg.second.plus
#
# data.2eq.det.408.0.gamma05 <- data.2eq
#
# saveRDS(gg,  paste0("objects/","gg.2eq.det.R250.p0.gamma05.rd"))
