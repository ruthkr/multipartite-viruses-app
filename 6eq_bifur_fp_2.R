library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions
R1.init <- 1
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

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 10000

get_steady_states <- function(num) {
	R1 <- numeric(num+1)
	gamma.h <- 1/num

	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							v1.init, v2.init,
							kappa1, kappa2, alpha, beta, 0, 0, sigma1, sigma2,
							epsilon1, epsilon2, delta1, delta2,
							num)

	options(scipen = 999)
	system("gcc -Ofast -lm 6eq_rk4_bipartite_ss.c -o compiled/rk4-ss")
	R1 <- system(paste("./compiled/rk4-ss", paste(params, collapse = ' '), "| cut -d, -f2"), intern = T)

	var.ss.fp2.df <- tibble(
		gamma = seq(0,1,gamma.h),
		R1 = as.numeric(R1)
	)
	return(var.ss.fp2.df)
}

var.ss.fp2.df <- get_steady_states(100000)

var.ss.fp2.df  %>%
	filter(R1 == 0) %>%
	head(1)

# Parameters --------------------------------------------------------------
gamma1 <- seq(0,1,0.00001)

a <- kappa1 * (sigma1/omega + 1)
b <- -(kappa1 + gamma1)
c <- gamma1

p.star2.plus <- (-b + sqrt(b^2 - 4*a*c)) / (2*a)
R1.star2.plus <- 1 - gamma1/(kappa1 * p.star2.plus)

p.star2.minus <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)
R1.star2.minus <- 1 - gamma1/(kappa1 * p.star2.minus)

ss.from.fp2 <- tibble(
	gamma1,
	p.plus = p.star2.plus,
	p.minus = p.star2.minus,
	R1.plus = R1.star2.plus,
	R1.minus = R1.star2.minus,
	R1.num = var.ss.fp2.df$R1
)



# ggplot(ss.from.fp2) +
# 	geom_line(aes(x = gamma1, y = R1.plus, colour = "R1.plus"), size=1) +
# 	geom_line(aes(x = gamma1, y = R1.minus, colour = "R1.minus"), linetype="dashed", size=1) +
# 	geom_line(aes(x = gamma1, y = R1.num, colour = "R1.num", group = (R1.num == 0)), linetype="dotted", size=1) +
# 	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
# 	scale_color_manual(values = c("R1.plus" = "red3", "R1.minus" = "royalblue","R1.num" = "orange" )) +
# 	theme_bw()

ggplot(ss.from.fp2) +
	geom_line(aes(x = gamma1, y = p.plus, colour = "p.plus"), size=1) +
	geom_line(aes(x = gamma1, y = p.minus, colour = "p.minus"), linetype="dashed", size=1) +
	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("p.plus" = "red3", "p.minus" = "royalblue")) +
	theme_bw()

