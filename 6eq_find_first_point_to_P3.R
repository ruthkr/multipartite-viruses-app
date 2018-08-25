library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 0 # Unused
R2.init <- 0 # Unused
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

gamma1 <- 0 # Unused
gamma2 <- gamma1 # Unused

sigma1 <- 0.1
sigma2 <- 0.1

epsilon1 <- 0.1
epsilon2 <- 0.1

delta1 <- 0.1
delta2 <- 0.1

# Partitions
gamma.part <- 55
R1.part <- 1000
R2.part <- 100

# Number of interations
options(scipen = 999)
n.iter = 25e3

# Concatenate parameters
params <- c(n.iter,
						R1.init, R2.init, p.init, s.init,
						v1.init, v2.init,
						kappa1, kappa2, alpha, beta,
						gamma1, gamma2, sigma1, sigma2,
						epsilon1, epsilon2, delta1, delta2,
						gamma.part, R1.part, R2.part)

# Compile
# system("gcc -Ofast -lm 6eq_rk4_bipartite_single_gamma_R1R2_tozero.c -o compiled/6eq-var-gamma")
# system("gcc -Ofast -lm 6eq_rk4_bipartite_var_gamma_R1R2_tozero.c -o compiled/6eq-var-gamma")
# system(paste("./compiled/6eq-var-gamma", paste(params, collapse = ' '), "> data/6eq_var_gamma_1000x100.csv"))

# Read simulation results -------------------------------------------------

# df <- data.table::fread(paste0("data/6eq_var_gamma.csv"))
df <- data.table::fread(paste0("data/6eq_var_gamma_100K.csv"))
# df.old <- data.table::fread(paste0("data/6eq_var_gamma.csv"))

clean_var_gamma_results_first_fp <- function(df) {
	df.clean <- df %>%
		dplyr::rename("gamma" = V1, "R1.init" = V2, "R2.init" = V3, "R1" = V4, "R2" = V5) %>%
		dplyr::filter(!(R1 == 0 & R2 == 0)) %>%
		group_by(gamma, R2.init) %>%
		summarise(R1.init = min(R1.init)) %>%
		select(gamma, R1.init, R2.init)

	return(df.clean)
}

clean_var_gamma_results_last_zero <- function(df) {
	df.clean <- df %>%
		dplyr::rename("gamma" = V1, "R1.init" = V2, "R2.init" = V3, "R1" = V4, "R2" = V5) %>%
		dplyr::filter((R1 == 0 & R2 == 0)) %>%
		group_by(gamma, R2.init) %>%
		summarise(R1.init = max(R1.init)) %>%
		select(gamma, R1.init, R2.init)

	return(df.clean)
}

df.clean.sep <- clean_var_gamma_results_last_zero(df)
df.clean.init <- clean_var_gamma_results_first_fp(df)

# df.clean.old <- clean_var_gamma_results_last_zero(df.old)
# df.clean.old <- df.clean.old[-284,]

# Compare results ------------------------------------------

compare <- F

if (compare) {
	identical(df.clean, df.clean.old)

	count <- 0
	for (i in 1:nrow(df.clean)) {
		if (!identical(df.clean.old[i,], df.clean[i,])) {
			print(paste("Rows", i, "are different."))
			print(paste("Old values:", toString(df.clean.old[i,])))
			print(paste("New values:", toString(df.clean[i,])))
			cat("\n")
			count <- count + 1
		}
	}

	df.clean.new <- data.frame(cbind(0,0,0)) %>%
		dplyr::rename("gamma" = X1, "R1.init" = X2, "R2.init" = X3)

	count <- 0
	for (i in 1:nrow(df.clean)) {
		if (!identical(df.clean.old[i,], df.clean[i,])) {
			count <- count + 1
			df.clean.new[count,] <- df.clean[i,]
		}
	}
}
