library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 1
R2.init <- 0.5
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

gamma1 <- 0.2
gamma2 <- 0.2

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 10000

# Flags
fixed.point.select <- "third.plus"

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
	gg2p <- gg +
		geom_hline( yintercept = p.star2.plus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star2.plus, color = "red3", linetype="dotted" )
} else if (fixed.point.select == "second.minus") {
	gg +
		geom_hline( yintercept = p.star2.minus, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star2.minus, color = "red3", linetype="dotted" )
}else if (fixed.point.select == "third.plus")  {
	gg +
		geom_hline( yintercept = p.star3.p, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star3.p, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star3.p, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star3.p, color = "black", linetype="dotted" ) +
		geom_hline( yintercept = v2.star3.p, color = "orange", linetype="dotted" )
}else if (fixed.point.select == "third.minus")  {
	gg +
		geom_hline( yintercept = p.star3.m, color = "springgreen4", linetype="dotted" ) +
		geom_hline( yintercept = R1.star3.m, color = "red3", linetype="dotted" ) +
		geom_hline( yintercept = R2.star3.m, color = "royalblue", linetype="dotted" ) +
		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" ) +
		geom_hline( yintercept = v1.star3.m, color = "black", linetype="dotted" ) +
		geom_hline( yintercept = v2.star3.m, color = "orange", linetype="dotted" )
}






# else if (fixed.point.select == "third.plus")  {
# 	gg +
# 		geom_hline( yintercept = p.star3.plus, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R1.star3.plus, color = "red3", linetype="dotted" ) +
# 		geom_hline( yintercept = R2.star3, color = "royalblue", linetype="dotted" ) +
# 		geom_hline( yintercept = s.star3.plus, color = "purple3", linetype="dotted" ) +
# 		geom_hline( yintercept = v1.star3.plus, color = "black", linetype="dotted" ) +
# 		geom_hline( yintercept = v2.star3.plus, color = "orange", linetype="dotted" )
# }else if (fixed.point.select == "third.minus")  {
# 	gg +
# 		geom_hline( yintercept = p.star3.minus, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R1.star3.minus, color = "red3", linetype="dotted" ) +
# 		geom_hline( yintercept = R2.star3, color = "royalblue", linetype="dotted" ) +
# 		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" ) +
# 		geom_hline( yintercept = v1.star3.minus, color = "black", linetype="dotted" ) +
# 		geom_hline( yintercept = v2.star3, color = "orange", linetype="dotted" )
# }else if (fixed.point.select == "third")  {
# 	gg +
# 		geom_hline( yintercept = p.star3, color = "springgreen4", linetype="dotted" ) +
# 		geom_hline( yintercept = R1.star3, color = "red3", linetype="dotted" ) +
# 		geom_hline( yintercept = R2.star3, color = "royalblue", linetype="dotted" ) +
# 		geom_hline( yintercept = s.star3, color = "purple3", linetype="dotted" ) +
# 		geom_hline( yintercept = v1.star3, color = "black", linetype="dotted" ) +
# 		geom_hline( yintercept = v2.star3, color = "orange", linetype="dotted" )
# }



# initial.values  <- tibble(R1.init = seq(0,1,0.2),
# 													R2.init = seq(0,1,0.2)) %>%
# 	expand(R1.init, R2.init)
#
# initial.values <- initial.values%>%
# 	mutate(value = numeric(nrow(initial.values)))
#
# alfR::lok_regar(for (combination in 1:nrow(initial.values)) {
#
# 	params <- c(n.iter,
# 							initial.values$R1.init[combination], initial.values$R2.init[combination],
# 							p.init, s.init, v1.init, v2.init,
# 							kappa1, kappa2, alpha, beta, gamma1, gamma2, sigma1, sigma2,
# 							epsilon1, epsilon2, delta1, delta2)
#
# 	system("gcc -Ofast -lm rk4_bipartite_6eq_test.c -o compiled/rk4-test")
# 	system(paste("./compiled/rk4-test", paste(params, collapse = ' '), " data/results-test.csv"))
#
# 	data <- data.table::fread("data/results-test.csv")
#
# 	initial.values$value[combination] = (data$V2 == data$V3) && (data$V2 == 0)
# })
#
# saveRDS(initial.values, "data/initial_values_matrix.RDS")
#
# initial.values <- initial.values %>%
# 	mutate( value = ifelse(value == 1, "to zero", "to fixed point" ))
#
#
# ggplot(initial.values) +
# 	geom_tile(aes(x = R1.init, y = R2.init, fill = as.factor(value))) +
# 	scale_fill_manual( values = c("to zero" = "mediumpurple1",
# 																"to fixed point" = "seagreen1") )

# library(tidyverse)
#
# initial.values <- readRDS("data/initial_values_matrix.RDS")
#
# # Function to find the slope and intercept -----------------
# find_line <- function(initial.values) {
# 	init.vals.cast <- reshape2::dcast(initial.values, R2.init ~ R1.init)
# 	rownames(init.vals.cast) <- init.vals.cast[,1]
# 	init.vals.cast[,1] <- NULL
#
# 	found.bot.left <- F
# 	found.top.right <- F
#
# 	dim <- nrow(init.vals.cast)
# 	parts <- dim - 1
#
# 	# R2 = 0 (bottom) x-axis
# 	if (!found.bot.left ) {
# 		for (j in 1:dim) {
# 			if (init.vals.cast[1,j] == "to fixed point") {
# 				bot.left <- c((j-1)/parts, 0)
# 				found.bot.left <- T
# 				break
# 			}
# 		}
# 	}
# 	# R1 = 0 (left) y-axis
# 	if (!found.bot.left ) {
# 		for (j in 1:dim) {
# 			if (init.vals.cast[j,1] == "to fixed point") {
# 				bot.left <- c(0, (j-1)/parts)
# 				found.bot.left <- T
# 				break
# 			}
# 		}
# 	}
#
# 	# R2 = 1 (top) x-axis
# 	if (!found.top.right ) {
# 		for (j in 1:dim) {
# 			if (init.vals.cast[dim,j] == "to fixed point") {
# 				top.right <- c((j-1)/parts, 1)
# 				found.top.right <- T
# 				break
# 			}
# 		}
# 	}
# 	# R1 = 1 (right) y-axis
# 	if (!found.top.right ) {
# 		for (j in 1:dim) {
# 			if (init.vals.cast[j,dim] == "to fixed point") {
# 				top.right <- c(1, (j-1)/parts)
# 				found.top.right <- T
# 				break
# 			}
# 		}
# 	}
# 	beta = (top.right[2] - bot.left[2]) / (top.right[1] - bot.left[1])
# 	alpha = bot.left[2] - beta * bot.left[1]
# 	return(c(alpha, beta))
# }
#
#
# # Run the function -----------------------------------------
#
# # Find the values
# line <- find_line(initial.values)
#
# ggplot(initial.values) +
# 	geom_tile(aes(x = R1.init, y = R2.init, fill = as.factor(value))) +
# 	scale_fill_manual( values = c("to zero" = "mediumpurple1",
# 																"to fixed point" = "seagreen1") ) +
# 	geom_abline(intercept = line[1], slope = line[2])
#
