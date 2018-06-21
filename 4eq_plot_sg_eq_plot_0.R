#for the case 4 eq and sigma1=sigma2=0

library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 1
R2.init <- 0.5
p.init <- 0.3
s.init <- 0.8

# Parameters
kappa1 <- 1
kappa2 <- 1
kappa <- 1

alpha <- 1
beta <- 1
omega <- 1

gamma1 <- 0.01
gamma2 <- 0.01
gamma <- 0.01

sigma1 <- 0
sigma2 <- 0

# Number of interations
n.iter <- 250

# Flags
num.equations <- 4
fixed.point.select <- "fourth"

# Put parameters in a vector
if (num.equations == 4) {
	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							kappa1, kappa2, alpha, beta, gamma1, gamma2, sigma1, sigma2)
} else if (num.equations == 6) {
	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							v1.init, v2.init,
							kappa1, kappa2, alpha, beta, gamma1, gamma2, sigma1, sigma2,
							epsilon1, epsilon2, delta1, delta2)
}

# Compile -----------------------------------------------------------------

if (num.equations == 4) {
	system("gcc -Ofast -lm rk4_bipartite_4eq_test.c -o compiled/rk4-test")
} else if (num.equations == 6) {
	system("gcc -Ofast -lm rk4_bipartite_6eq_test.c -o compiled/rk4-test")
}

system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

# Read results and plot ---------------------------------------------------

# Read data
data <- data.table::fread("data/results-test.csv")

# Create plot
gg <- ggplot(data) +
	geom_line(aes(x = V1, y = V2, colour = "R1")) +
	geom_line(aes(x = V1, y = V3, colour = "R2")) +
	geom_line(aes(x = V1, y = V4, colour = "p")) +
	geom_line(aes(x = V1, y = V5, colour = "S")) +
	labs(x = "Step", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R1" = "red3",        "R2" = "royalblue",
																"p" = "springgreen4", "S" = "purple3")) +
	theme_bw()


# First Fixed Point
p.star1 <- p.init
s.star1 <- s.init

# Second Fixed Point
p.star2 <- tail(data$V4,1)
R2.star2 <- 1 - gamma2/ (kappa2*p.star2)
s.star2 <- 1 - p.star2

# Third Fixed Point
s.star3 <- s.init
p.star3 <- 1 - s.star3
R1.star3 <- 1 - gamma1/ (kappa1*p.star3)


# Fourth Fixed Point
p.star4 <- tail(data$V4,1)
R1.star4 <- tail(data$V2,1)
R2.star4 <- 1 - R1.star4 - gamma2 / (kappa2*p.star4)
s.star4 <- 1 - p.star4


if (fixed.point.select == "first") {
	gg +
		geom_hline( yintercept = p.star1, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = s.star1, color = "purple3", linetype="dotted" )
} else if (fixed.point.select == "second") {
	gg +
		geom_hline( yintercept = p.star2, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R2.star2, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star2, color = "purple3", linetype="dotted" )
}else if (fixed.point.select == "third") {
	gg +
		geom_hline( yintercept = p.star3, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star3, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" )
}else if (fixed.point.select == "fourth")  {
	gg +
		geom_hline( yintercept = p.star4, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star4, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star4, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star4, color = "purple3", linetype="dotted" )
}

