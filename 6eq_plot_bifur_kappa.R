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
alpha <- 1
beta <- 1
omega <- 1

gamma1 <- 0.0112866
gamma2 <- 0.0112866

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 7500

# Flags
fixed.point.select <- "third"

#Compile
compile <- function(num, point = "third") {
	kappa.h <- 1/num
	R1.ss <- numeric(num)
	R2.ss <- numeric(num)
	i <- 1

	if (point == "second") {
		s.init <- 0
		R2.init <- 0
	}

	for (kappa in seq(0,1,kappa.h)) {
		params <- c(n.iter,
								R1.init, R2.init, p.init, s.init,
								v1.init, v2.init,
								kappa, kappa, alpha, beta, gamma1, gamma2, sigma1, sigma2,
								epsilon1, epsilon2, delta1, delta2)

		system("gcc -Ofast -lm rk4_bipartite_6eq_test.c -o compiled/rk4-test")
		system(paste("./compiled/rk4-test", paste(params, collapse = ' '), "> data/results-test.csv"))

		# Read data
		data <- data.table::fread("data/results-test.csv")

		# Store steady state value of R1
		R1.ss[i] <- data$V2[nrow(data)]
		R2.ss[i] <- data$V3[nrow(data)]
		i = i + 1
	}
	return(tibble(R1.ss, R2.ss))
}

get_steady_states <- function(num) {
	R1.3.ss <- compile(num)$R1.ss
	R2.3.ss <- compile(num)$R2.ss
	R1.2.ss <- compile(num, "second")$R1.ss
	R2.2.ss <- compile(num, "second")$R2.ss

	kappa.h <- 1/num
	RNA.ss.df <- tibble(
		kappa = seq(0,1,kappa.h),
		R1.3 = R1.3.ss,
		R1.2 = R1.2.ss,
		R2.3 = R2.3.ss,
		R2.2 = R2.2.ss
	)
	return(RNA.ss.df)
}

alfR::lok_regar( RNA.ss.df <- get_steady_states(1000) )

RNA.ss.df %>%
	filter(R1.2 == 0) %>%
	head(1)

# Create plot
ggplot(RNA.ss.df) +
	geom_line(aes(x = kappa, y = R1.2, colour = "R1.3")) +
	geom_line(aes(x = kappa, y = R2.2, colour = "R2.3")) +
	labs(x = "Kappa", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("R1.3" = "red3", "R2.3" = "mediumpurple3")) +
	theme_bw()


RNA.data.slice <- RNA.ss.df %>%
	slice(50:60) %>%
	select(gamma, R1.2, R2.2)
