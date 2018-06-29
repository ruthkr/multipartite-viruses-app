library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 0.75
R2.init <- 0
p.init <- 0
s.init <- 0
v1.init <- 0
v2.init <- 0

# Parameters
kappa1 <- 1
kappa2 <- 1

alpha <- 1
beta <- 1
omega <- 1

gamma1 <- 0.45
gamma2 <- 0.45

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 10000

# Flags
fixed.point.select <- "second.plus"

params <- c(n.iter,
						R1.init, R2.init, p.init, s.init,
						v1.init, v2.init,
						kappa1, kappa2, alpha, beta, gamma1, gamma2, sigma1, sigma2,
						epsilon1, epsilon2, delta1, delta2)

#Compile
system("gcc -Ofast -lm 6eq_rk4_bipartite_test.c -o compiled/rk4-test")
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
	labs(x = "Time", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("R1" = "red3",        "R2" = "royalblue",
																"p" = "springgreen4", "S" = "purple3",
																"v1" = "black", "v2" = "orange")) +
	theme_bw()

tail(data)


# Second Fixed Point
a <- kappa1 * (sigma1/omega + 1)
b <- -(kappa1 + gamma1)
c <- gamma1

p.star2.plus <- (-b + sqrt(b^2 - 4*a*c)) / (2*a)
R1.star2.plus <- 1 - gamma1/(kappa1 * p.star2.plus)

p.star2.minus <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)
R1.star2.minus <- 1 - gamma1/(kappa1 * p.star2.minus)

# Third Fixed Point
s.star3 <- tail(data$V5,1)

a3 <- (sigma1 / omega + 1 ) * kappa1
b3 <- (sigma1 / omega) * kappa1 * s.star3 - kappa1 * (1 - s.star3) - epsilon1 * s.star3 - gamma1
c3 <- (epsilon1*s.star3 + gamma1) * (1 - s.star3)

p.star3.p <- (-b3 + sqrt(b3^2 - 4*a3*c3)) / (2*a3)
R1.star3.p <- (kappa1*p.star3.p - epsilon1*s.star3 - gamma1) / (kappa1 * (p.star3.p + s.star3))
R2.star3.p <- R1.star3.p * s.star3 / p.star3.p
v1.star3.p <- epsilon1*R1.star3.p*s.star3 / delta1
v2.star3.p <- epsilon2*R2.star3.p*s.star3 / delta2

p.star3.m <- (-b3 - sqrt(b3^2 - 4*a3*c3)) / (2*a3)
R1.star3.m <- (kappa1*p.star3.m - epsilon1*s.star3 - gamma1) / (kappa1 * (p.star3.m + s.star3))
R2.star3.m <- R1.star3.m * p.star3.m / s.star3
v1.star3.m <- epsilon1*R1.star3.m*s.star3 / delta1
v2.star3.m <- epsilon2*R2.star3.m*s.star3 / delta2

if (fixed.point.select == "second.plus") {
	gg +
		geom_hline( yintercept = p.star2.plus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star2.plus, color = "red3", linetype="dotted" )
} else if (fixed.point.select == "second.minus") {
	gg +
		geom_hline( yintercept = p.star2.minus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star2.minus, color = "red3", linetype="dotted" )
}else if (fixed.point.select == "third.plus")  {
	gg +
		geom_hline( yintercept = p.star3.plus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star3.plus, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star3, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star3.plus, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star3.plus, color = "black", linetype="dotted" ) +
		geom_hline( yintercept = v2.star3.plus, color = "orange", linetype="dotted" )
}else if (fixed.point.select == "third.minus")  {
	gg +
		geom_hline( yintercept = p.star3.m, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star3.m, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star3.m, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star3.m, color = "black", linetype="dotted" ) +
		geom_hline( yintercept = v2.star3.m, color = "orange", linetype="dotted" )
}




