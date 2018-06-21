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

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Number of interations
n.iter <- 10000

# Flags
fixed.point.select <- "second"

#Compile
compile <- function(num, point = "third") {
	gamma.h <- 1/num
	R1.ss <- numeric(num)
	R2.ss <- numeric(num)
	i <- 1

	if (point == "second") {
		s.init <- 0
		R2.init <- 0
	}

	for (gamma in seq(0,1,gamma.h)) {
		params <- c(n.iter,
								R1.init, R2.init, p.init, s.init,
								v1.init, v2.init,
								kappa1, kappa2, alpha, beta, gamma, gamma, sigma1, sigma2,
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
	# R1.2.ss <- compile(num, "second")$R1.ss
	# R2.2.ss <- compile(num, "second")$R2.ss


	gamma.h <- 1/num
	RNA.ss.df <- tibble(
		gamma = seq(0,1,gamma.h),
		R1.3 = R1.3.ss,
		# R1.2 = R1.2.ss,
		R2.3 = R2.3.ss
		# R2.2 = R2.2.ss
		)
	return(RNA.ss.df)
}

alfR::lok_regar( RNA.ss.df <- get_steady_states(100) )

RNA.ss.df %>%
	filter(R1.3 == R2.3) %>%
	filter(R1.3 == 0) %>%
	head(1)

# Create plot
gR1 <- ggplot(RNA.ss.df) +
	geom_line(aes(x = gamma, y = R1.3, colour = "R1.3")) +
	geom_line(aes(x = gamma, y = R1.2, colour = "R1.2")) +
	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("R1.3" = "red3", "R1.2" = "royalblue")) +
	theme_bw()

gR2 <- ggplot(RNA.ss.df) +
	geom_line(aes(x = gamma, y = R2.3, colour = "R2.3")) +
	geom_line(aes(x = gamma, y = R2.2, colour = "R2.2")) +
	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("R2.3" = "mediumpurple3", "R2.2" = "lightpink2")) +
	theme_bw()

gR1.R2.fp.3 <- ggplot(RNA.ss.df) +
	geom_line(aes(x = gamma, y = R1.3, colour = "R1.3")) +
	geom_line(aes(x = gamma, y = R2.3, colour = "R2.3")) +
	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("R1.3" = "red3", "R2.3" = "royalblue")) +
	theme_bw()

gR1.R2.fp.2 <- ggplot(RNA.ss.df) +
	geom_line(aes(x = gamma, y = R2.2, colour = "R2.2")) +
	geom_line(aes(x = gamma, y = R1.2, colour = "R1.2")) +
	labs(x = "Gamma", y = "Equilibrium Point", color = "Variable") +
	scale_color_manual(values = c("R1.2" = "red3", "R2.2" = "royalblue")) +
	theme_bw()

# gR1.R2.fp.2 + ggsave(filename = "bifur_fp2_kappa1.pdf", width = 5.75, height = 3, dpi = 96, device = cairo_pdf)

# gR1.R2.fp.3 + ggsave(filename = "bifur_fp3_kappa1.pdf", width = 5.75, height = 3, dpi = 96, device = cairo_pdf)
