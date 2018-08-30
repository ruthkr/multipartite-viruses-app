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

epsilon1 <- 0
epsilon2 <- 0

delta1 <- 0
delta2 <- 0

# Number of interations
n.iter <- 10000

get_steady_states <- function(num) {
	R1 <- numeric(num+1)
	p <- numeric(num+1)
	gamma.h <- 1/num

	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							v1.init, v2.init,
							kappa1, kappa2, alpha, beta, 0, 0, sigma1, sigma2,
							epsilon1, epsilon2, delta1, delta2,
							num)

	options(scipen = 999)
	system("gcc -Ofast -lm 6eq_rk4_bipartite_ss.c -o compiled/rk4-ss")
	temp <- system(paste("./compiled/rk4-ss", paste(params, collapse = ' '), "| cut -d, -f2,4"), intern = T)
	temp <- str_split_fixed(temp, ", ", 2)

	var.ss.fp2.df <- tibble(
		gamma = seq(0,1,gamma.h),
		R1 = as.numeric(temp[,1]),
		p = as.numeric(temp[,2])
		)
	return(var.ss.fp2.df)
}

var.ss.fp2.df <- get_steady_states(1000)

var.ss.fp2.df  %>%
	filter(R1 == 0) %>%
	head(1)

# Parameters --------------------------------------------------------------
gamma1 <- seq(0, 1, 0.001)

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
	p.num = var.ss.fp2.df$p,
	R1.plus = R1.star2.plus,
	R1.minus = R1.star2.minus,
	R1.num = var.ss.fp2.df$R1
)

# gg <- ggplot(ss.from.fp2) +
# 	geom_line(aes(x = gamma1, y = R1.plus, colour = "R1.plus"), size=1) +
# 	geom_line(aes(x = gamma1, y = R1.minus, colour = "R1.minus"), linetype="dashed", size=1) +
# 	geom_line(aes(x = gamma1, y = R1.num, colour = "R1.num", group = (R1.num == 0)), linetype="dotdash", size=1) +
# 	labs(x = "Gamma", y = "Equilibrium point", color = "Variable") +
# 	# annotate("text", x = 0.15, y = 0.5, label = "Bistable\n system") +
# 	# annotate("text", x = 0.80, y = 0.5, label = "Monostable\n system") +
# 	# Lower arrows
# 	geom_segment(aes(x = 0.65, y = 0.18, xend = 0.65, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.85, y = 0.18, xend = 0.85, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.4, y = 0.07, xend = 0.4, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.5, y = 0.12, xend = 0.5, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	# Upper arrows
# 	geom_segment(aes(x = 0.2, y = 0.99, xend = 0.2, yend = 0.83), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.35, y = 0.81, xend = 0.35, yend = 0.65), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.5, y = 0.59, xend = 0.5, yend = 0.43), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	# Inner arrows up
# 	geom_segment(aes(x = 0.05, y = 0.72, xend = 0.05, yend = 0.88), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.25, y = 0.49, xend = 0.25, yend = 0.65), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	# Inner arrows down
# 	geom_segment(aes(x = 0.05, y = 0.06, xend = 0.05, yend = 0.22), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.25, y = 0.09, xend = 0.25, yend = 0.25), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	geom_segment(aes(x = 0.45, y = 0.2, xend = 0.45, yend = 0.36), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
# 	scale_color_manual(breaks = c("R1.plus", "R1.minus","R1.num"),
# 										 values = c("R1.plus" = "red3",
# 																"R1.minus" = "royalblue",
# 																"R1.num" = "darkgoldenrod1"),
# 										 labels = c("R1.plus" = expression(paste(R["1+     "])),
# 										 					 "R1.minus" = expression(paste(R["1-      "])),
# 										 					 "R1.num" = expression(paste(R["1num"])))) +
# 	scale_x_continuous(breaks = seq(0, 1, 0.1)) +
# 	theme_bw()
#
# gg
#
# saveRDS(gg, "objects/plot_bifur_2eq_without_annotate.rds")



ggp <- ggplot(ss.from.fp2) +
	geom_line(aes(x = gamma1, y = p.plus, colour = "p.plus"), size=1) +
	geom_line(aes(x = gamma1, y = p.minus, colour = "p.minus"), linetype="dashed", size=1) +
	geom_line(aes(x = gamma1, y = p.num, colour = "p.num", group = (p.num == 0)), linetype="dotdash", size=1) +
	labs(x = "Gamma", y = "Fixed Point", color = "Variable") +
	# Lower arrows
	# geom_segment(aes(x = 0.7, y = 0.18, xend = 0.7, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.85, y = 0.18, xend = 0.85, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 1, y = 0.18, xend = 1, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.55, y = 0.18, xend = 0.55, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.4, y = 0.18, xend = 0.4, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.25, y = 0.18, xend = 0.25, yend = 0.02), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.5, y = 0.55, xend = 0.5, yend = 0.44), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.375, y = 0.37, xend = 0.375, yend = 0.26), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# Upper arrows
	# geom_segment(aes(x = 0.2, y = 1, xend = 0.2, yend = 0.91), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.35, y = 0.96, xend = 0.35, yend = 0.87), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.5, y = 0.91, xend = 0.5, yend = 0.8), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# Inner arrows up
	# geom_segment(aes(x = 0.05, y = 0.72, xend = 0.05, yend = 0.88), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.25, y = 0.69, xend = 0.25, yend = 0.85), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.45, y = 0.62, xend = 0.45, yend = 0.78), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# Inner arrows down
	# geom_segment(aes(x = 0.05, y = 0.12, xend = 0.05, yend = 0.28), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
	# geom_segment(aes(x = 0.25, y = 0.32, xend = 0.25, yend = 0.48), size = 0.4, arrow = arrow(length = unit(0.2, "cm"))) +
		scale_color_manual(breaks = c("p.plus", "p.minus","p.num"),
											 values = c("p.plus" = "red3",
																	"p.minus" = "royalblue",
																	"p.num" = "darkgoldenrod1"),
											 labels = c("p.plus" = expression(paste(p["+     "])),
											 					 "p.minus" = expression(paste(p["-      "])),
											 					 "p.num" = expression(paste(p["num"])))) +
		scale_x_continuous(breaks = seq(0, 1, 0.1)) +
	theme_bw()

saveRDS(ggp, "objects/plot_bifur_2eq_p_without_annotate.rds")
