library(tidyverse)

# Parameters --------------------------------------------------------------
# Initial conditions

R1.init <- 0 # Unused
p.init <- 0 # Unused

# Parameters
kappa <- 1
omega <- 1
gamma <- 0 # Unused
sigma <- 0.1

# Partitions
gamma.part <- 55
R1.part <- 1000
p.part <- 10

# Number of interations
options(scipen = 999)
n.iter = 50000

# Concatenate parameters
params <- c(n.iter,
						R1.init, p.init,
						kappa, omega, gamma, sigma,
						gamma.part, R1.part, p.part)

# Compile
# system("gcc -Ofast -lm 2eq_rk4_rep_var_gamma_tozero.c -o compiled/2eq-var-gamma")
# system(paste("./compiled/2eq-var-gamma", paste(params, collapse = ' '), "> data/2eq_var_gamma.csv"))

# Read simulation results -------------------------------------------------

df <- data.table::fread(paste0("data/2eq_var_gamma.csv"))
# df.old <- data.table::fread(paste0("data/2eq_var_gamma_old.csv"))

clean_var_gamma_results <- function(df) {
	df.clean <- df %>%
		dplyr::rename("gamma" = V1, "R1.init" = V2, "p.init" = V3, "R1" = V4, "p" = V5) %>%
		dplyr::filter(!(R1 == 0 & p == 0)) %>%
		group_by(gamma, p.init) %>%
		summarise(R1.init = min(R1.init)) %>%
		select(gamma, R1.init, p.init)

	return(df.clean)
}

df.clean <- clean_var_gamma_results(df)
# df.clean.old <- clean_var_gamma_results(df.old)


# Compare results ------------------------------------------

compare <- FALSE

if (compare) {
	identical(df.clean, df.clean.old)

	for (i in 1:nrow(df.clean)) {
		if (!identical(df.clean.old[i,], df.clean[i,])) {
			print(paste("Rows", i, "are different."))
			print(paste("Old values:", toString(df.clean.old[i,])))
			print(paste("New values:", toString(df.clean[i,])))
			cat("\n")
		}
	}

	df.clean.new <- data.frame(cbind(0,0,0)) %>%
		dplyr::rename("gamma" = X1, "R1.init" = X2, "p.init" = X3)

	count <- 0
	for (i in 1:nrow(df.clean)) {
		if (!identical(df.clean.old[i,], df.clean[i,])) {
			count <- count + 1
			df.clean.new[count,] <- df.clean[i,]
		}
	}
}
