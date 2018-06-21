library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 0.038
R2.init <- 0.02
p.init <- 0.1
s.init <- 0.1


# Parameters
kappa1 <- 0.92
kappa2 <- 0.92


alpha <- 1
beta <- 1
omega <- 1

gamma1 <- 0.75
gamma2 <- 0.75
gamma <- 0.75

sigma1 <- 0.01
sigma2 <- 0.01


# Number of interations
n.iter <- 1000

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
																"p" = "springgreen4", "S" = "purple3",
																"v1" = "black", "v2" = "orange")) +
	theme_bw()

#Second Fixed Point
R2.star2 <- tail(data$V3,1)
s.star2 <- omega*R2.star2 / (omega*R2.star2 + sigma2)

# Third fixed point
p.star.plus <- ( (kappa1 + gamma1) + sqrt( (kappa1 + gamma1)^2 - 4 * kappa1 * gamma1 * ( (sigma1 / alpha) + 1) ) ) / ( 2 * kappa1 * ( (sigma1 / alpha) + 1) )
p.star.minus <- ( (kappa1 + gamma1) - sqrt( (kappa1 + gamma1)^2 - 4 * kappa1 * gamma1 * ( (sigma1 / alpha) + 1) ) ) / ( 2 * kappa1 * ( (sigma1 / alpha) + 1) )

R1.star.plus <- 1 - gamma1 / (kappa1 * p.star.plus)
R1.star.minus <- 1 - gamma1 / (kappa1 * p.star.minus)

# Fourth Fixed Points
# R1.star <- tail(data,1)$V2
# p.star <- (omega*R1.star*sigma2 + omega*sigma1*(gamma/kappa)) / (sigma1*sigma2 + omega*sigma1*(1-R1.star) + omega*R1.star*sigma2)
# R2.star <- (1 - R1.star) - (gamma / (kappa*p.star))
# s.star <- omega*((1-R1.star) - gamma/(kappa*p.star)) * (1 - p.star) / (sigma2 + omega * ((1 - R1.star) - gamma / (kappa * p.star)))

R1.star <- tail(data,1)$V2
p.star <- (omega*R1.star*sigma2 + omega*sigma1*(gamma2/kappa2)) / (sigma1*sigma2 + omega*sigma1*(1-R1.star) + omega*R1.star*sigma2)
R2.star <- (1 - R1.star) - (gamma2 / (kappa2*p.star))
s.star <- omega*((1-R1.star) - gamma2/(kappa2*p.star)) * (1 - p.star) / (sigma2 + omega * ((1 - R1.star) - gamma2 / (kappa2 * p.star)))


if (fixed.point.select == "second") {
	gg +
		geom_hline( yintercept = s.star2, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star2, color = "royalblue", linetype="dotted" )
} else if (fixed.point.select == "plus.third") {
	gg +
		geom_hline( yintercept = p.star.plus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star.plus, color = "red3", linetype="dotted" )
} else if (fixed.point.select == "minus.third") {
	gg +
		geom_hline( yintercept = p.star.minus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star.minus, color = "red3", linetype="dotted" )
}else if (fixed.point.select == "fourth")  {
	gg +
		geom_hline( yintercept = p.star, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star, color = "purple3", linetype="dotted" )
}


# R2.star <- tail(data,1)$V3
# a <- kappa1 * ((sigma1/ alpha) * (sigma2 + beta * R2.star) + (1 - R2.star) * sigma2)
# b <- - sigma2 * ( kappa1 * ( 1 - R2.star ) + gamma1)
# c <- gamma1 * sigma2
# p.star.plus4 <- (- b + sqrt ( b^2 - 4*a*c)) / 2*a
# p.star.minus4 <- (- b - sqrt ( b^2 - 4*a*c)) / 2*a
#
# R1.star.plus4 <- (1 - R2.star) - (gamma1 / (kappa1*p.star.plus4))
# R1.star.minus4 <- (1 - R2.star) - (gamma1 / (kappa1*p.star.minus4))
#
# s.star.plus4 <- beta * R2.star * (1 - p.star.plus4) / (beta*R2.star + sigma2)
# s.star.minus4 <- beta * R2.star * (1 - p.star.minus4) / (beta*R2.star + sigma2)
#
#
# if (fixed.point.select == "plus.fourth") {
# 	gg +
# 		geom_hline( yintercept = p.star.plus4, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R1.star.plus4, color = "red3", linetype="dotted" ) +
# 		geom_hline( yintercept = R2.star, color = "royalblue", linetype="dotted" ) +
# 		geom_hline( yintercept = s.star.plus4, color = "purple3", linetype="dotted" )
# } else if (fixed.point.select == "minus.fourth") {
# 	gg +
# 		geom_hline( yintercept = p.star.minus4, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R1.star.minus4, color = "red3", linetype="dotted" ) +
# 		geom_hline( yintercept = R2.star, color = "royalblue", linetype="dotted" ) +
# 		geom_hline( yintercept = s.star.minus4, color = "purple3", linetype="dotted" )
# }
