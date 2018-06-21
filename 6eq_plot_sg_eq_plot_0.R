library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 1
R2.init <- 1
p.init <- 0
s.init <- 0
v1.init <- 0
v2.init <- 0

# Parameters
kappa1 <- 1
kappa2 <- 1
kappa <- 1

alpha <- 1
beta <- 1
omega <- 1

gamma1 <- 0.9
gamma2 <- 0.9
gamma <- 0.9

sigma1 <- 0
sigma2 <- 0

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 10000

# Flags
fixed.point.select <- "third"


params <- c(n.iter,
						R1.init, R2.init, p.init, s.init,
						v1.init, v2.init,
						kappa1, kappa2, alpha, beta, gamma1, gamma2, sigma1, sigma2,
						epsilon1, epsilon2, delta1, delta2)

#Compile
system("gcc -Ofast -lm rk4_bipartite_6eq_test.c -o compiled/rk4-test")
system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

# Read data
data <- data.table::fread("data/results-test.csv")

# Create plot
gg <- ggplot(data) +
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "R2")) +
	geom_line(aes(x = V1, y = V4, colour = "p")) +
	geom_line(aes(x = V1, y = V5, colour = "S")) +
	geom_line(aes(x = V1, y = V6, colour = "v1")) +
	geom_line(aes(x = V1, y = V7, colour = "v2")) +
	labs(x = "Step", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R1" = "red3",        "R2" = "royalblue",
																"p" = "springgreen4", "S" = "purple3",
																"v1" = "black", "v2" = "orange")) +
	theme_bw()

# First Fixed Point
s.star1 <- tail(data$V5,1)
p.star1 <- 1 - s.star1
R1.star1 <- 1 - epsilon1*s.star1 / (kappa1*p.star1) - gamma1 / (kappa1*p.star1)
v1.star1 <- epsilon1*R1.star1*s.star1 / delta1

# Second Fixed Point
s.star2 <- tail(data$V5,1)
p.star2 <- 1 - s.star2
R2.star2 <- 1 - epsilon2*s.star2 / (kappa2*p.star2) - gamma2 / (kappa2*p.star2)
v2.star2 <- epsilon2*R2.star2*s.star2 / delta2

# Third Fixed Point
R1.star3 <- tail(data$V2,1)
p.star3 <- tail(data$V4,1)
s.star3 <- tail(data$V5,1)
R2.star3 <- 1 - R1.star3 - epsilon2*s.star3 / (kappa2*p.star3) - gamma2 / (kappa2*p.star3)
v1.star3 <- epsilon1*R1.star3*s.star3 / delta1
v2.star3 <- epsilon2*R2.star3*s.star3 / delta2

# Four Fixed Point
R1.star4 <- R1.init
p.star4 <- p.init
s.star4 <- s.init
R2.star4 <- R2.init
v1.star4 <- v1.init
v2.star4 <- v2.init

if (fixed.point.select == "second") {
	gg +
		geom_hline( yintercept = p.star2, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R2.star2, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star2, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v2.star2, color = "orange", linetype="dotted" )
} else if (fixed.point.select == "first") {
	gg +
		geom_hline( yintercept = p.star1, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star1, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star1, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star1, color = "orange", linetype="dotted" )
}else if (fixed.point.select == "third")  {
	gg +
		geom_hline( yintercept = p.star3, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star3, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star3, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star3, color = "black", linetype="dotted" ) +
		geom_hline( yintercept = v2.star3, color = "orange", linetype="dotted" )
}else if (fixed.point.select == "fourth")  {
	gg +
		geom_hline( yintercept = p.star4, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star4, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star4, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star4, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star4, color = "black", linetype="dotted" ) +
		geom_hline( yintercept = v2.star4, color = "orange", linetype="dotted" )
}
