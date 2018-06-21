library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R.init <- 1
p.init <- 0

# Parameters
kappa <- 1
omega <- 1
gamma <- 0.1
sigma <- 0.1

# Number of interations
n.iter <- 1000

# Flags
fixed.point.select <- "second.plus"

params <- c(n.iter,
						R.init, p.init,
						kappa, omega, gamma, sigma)

#Compile
system("gcc -Ofast -lm rk4_bipartite_2eq_test.c -o compiled/rk4-test")
system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

# Read data
data.2eq <- data.table::fread("data/results-test.csv")

# Create plot
gg <- ggplot(data.2eq) +
	geom_line(aes(x = V1, y = V2, colour = "RNA")) +
	geom_line(aes(x = V1, y = V3, colour = "viral.replicase")) +
	labs(x = "Time", y = "Variables", color = "Variable") +
	scale_color_manual(values = c("RNA" = "red3",
																"viral.replicase" = "springgreen4")) +
	theme_bw()


# Second Fixed Point
a <- kappa * (sigma/omega + 1)
b <- -(kappa + gamma)
c <- gamma

p.plus <- (-b + sqrt(b^2 - 4*a*c)) / (2*a)
R.plus <- 1 - gamma/(kappa * p.plus)

p.minus <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)
R.minus <- 1 - gamma/(kappa * p.minus)


if (fixed.point.select == "second.plus") {
  gg +
		geom_hline( yintercept = p.plus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R.plus, color = "red3", linetype="dotted" )
} else if (fixed.point.select == "second.minus") {
	gg +
		geom_hline( yintercept = p.minus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R.minus, color = "red3", linetype="dotted" )
}



