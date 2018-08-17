# Bifurcation Diagram obtained from the fixed points

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

alpha <- 1
beta <- 1

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 10000

get_steady_states <- function(num) {
	R2 <- numeric(num+1)
	s <- numeric(num+1)
	gamma.h <- 1/num

	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							v1.init, v2.init,
							kappa1, kappa2, alpha, beta, 0, 0, sigma1, sigma2,
							epsilon1, epsilon2, delta1, delta2,
							num)

	options(scipen = 999)
	system("gcc -Ofast -lm 6eq_rk4_bipartite_ss.c -o compiled/rk4-ss")
	temp <- system(paste("./compiled/rk4-ss", paste(params, collapse = ' '), "| cut -d, -f2,3,5"), intern = T)
	temp <- str_split_fixed(temp, ", ", 3)

	var.ss.fp3.df <- tibble(
		gamma = seq(0,1,gamma.h),
		R1 = as.numeric(temp[,1]),
		R2 = as.numeric(temp[,2]),
		s = as.numeric(temp[,3])
	)
	return(var.ss.fp3.df)
}

var.ss.fp3.df <- get_steady_states(10000)

# var.ss.fp3.df  %>%
# 	filter(R1 == 0) %>%
# 	head(1)


# Plot the bifurcation manually using the equation of fixed point
omega <- 1
gamma <- var.ss.fp3.df$gamma


# var.ss.fp3.df.bis <- var.ss.fp3.df %>%
# 	mutate(
# 		a = (sigma1 / omega + 1 ) * kappa1,
# 		b = (sigma1 / omega) * kappa1 * s - kappa1 * (1 - s) - epsilon1 * s - gamma,
# 		c = (epsilon1*s + gamma) * (1 - s),
# 		p.star.plus = (-b + sqrt(b^2 - 4 * a * c)) / (2*a),
# 		p.star.minus = (-b - sqrt(b^2 - 4 * a * c)) / (2*a),
# 		R1.star.plus = (kappa1*p.star.plus - epsilon1*s - gamma) / (kappa1 * (p.star.plus + s)),
# 		R1.star.minus = (kappa1*p.star.minus - epsilon1*s - gamma) / (kappa1 * (p.star.minus + s))
# 	) %>%
# 	select(-a,-b,-c)


var.ss.fp3.df.bis <- var.ss.fp3.df %>%
	mutate(
		a = (sigma1 / omega + 1 ) * kappa1,
		b = (sigma1 / omega) * kappa1 * s - kappa1 * (1 - s) - epsilon1 * s - gamma,
		c = (epsilon1*s + gamma) * (1 - s),
		p.star.plus = (-b + sqrt(b^2 - 4 * a * c)) / (2*a),
		p.star.minus = (-b - sqrt(b^2 - 4 * a * c)) / (2*a),
		R2.star.plus = sigma2 * s / (omega * (1 - p.star.plus - s)),
		R1.star.plus = 1 - R2.star.plus - epsilon1*s / (kappa1*p.star.plus) - gamma / (kappa1*p.star.plus),
		R2.star.minus = sigma2 * s / (omega * (1 - p.star.minus - s)),
		R1.star.minus = 1 - R2.star.minus - epsilon1*s / (kappa1*p.star.minus) - gamma / (kappa1*p.star.minus)
	) %>%
	select(-a,-b,-c)



# Create plot

# ggplot(var.ss.fp3.df.bis) +
# 	 	geom_line(aes(x = gamma, y = p.star.plus, colour = "p.plus"), size=1) +
# 	 	geom_line(aes(x = gamma, y = p.star.minus, colour = "p.minus"), linetype="dashed", size=1) +
# 	 	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
# 	 	scale_color_manual(values = c("p.plus" = "red3", "p.minus" = "royalblue")) +
# 	 	theme_bw()

ggplot(var.ss.fp3.df.bis) +
	geom_line(aes(x = gamma, y = R1.star.plus, colour = "R1.plus"), size=1) +
	geom_line(aes(x = gamma, y = R1.star.minus, colour = "R1.minus"), linetype="dashed", size=1) +
	# geom_line(aes(x = gamma, y = R1, colour = "R1"), linetype="dashed", size=1) +
	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("R1.plus" = "red3", "R1.minus" = "royalblue", "R1" = "yellow")) +
	scale_x_continuous(breaks = seq(0, 1, 0.05)) +
	theme_bw()


var.ss.fp3.df.bis  %>%
	filter(s == 0) %>%
	head(1)
