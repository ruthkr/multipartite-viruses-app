library(tidyverse)

source("6eq_find_first_point_to_P3.R")

# Compile
# system("gcc -Wall -Ofast -lm 6eq_ssa_prob_exit.c -o compiled/6eq_prob_exit")

# General parameters
max.time <- 1000
n.sim <- 1000

df.clean <- df.clean.init

scale <- 1000
max.run <- nrow(df.clean)
# first.run <- 1
first.run <- 2882
len.run <- max.run - first.run + 1

curr.capacity <- 1000


# Run simulations (standard) ----------------------------------------------

# Create data frame
# results <- data.frame(cbind(0,0,0,0))
# results <- results %>%
# 	dplyr::rename("gamma" = X1, "R1.init" = X2, "R2.init" = X3, "prob" = X4)
#
# alfR::lok_regar(
# 	for (i in 1:n.run) {
# 		print(paste0("Simulation ", i, "/", n.run, " started at ", format(Sys.time(), "%H:%M:%S")))
# 		# Read parameters from table
# 		gamma <- as.numeric(df.clean[i, "gamma"])
# 		R1.init <- as.integer(as.numeric(df.clean[i, "R1.init"]) * scale)
# 		R2.init <- as.integer(as.numeric(df.clean[i, "R2.init"]) * scale)
#
# 		# Update parameters
# 		params <- c(max.time,
# 								R1.init, R2.init, 0, 0, 0, 0,
# 								kappa1, omega, gamma, sigma1,
# 								epsilon1, delta1, curr.capacity,
# 								n.sim)
# 		# Run simulation
# 		prob <- as.numeric(system(paste("./compiled/6eq_prob_exit", paste(params, collapse = ' ')), intern = T))
# 		results[1,] <- cbind(gamma, R1.init, R2.init, prob)
#
# 		# print(results)
# 		if (i == 1) {
# 			write_csv(results, "data/6eq_ssa_prob_extinction_old.csv")
# 		} else {
# 			write_csv(results, "data/6eq_ssa_prob_extinction_old.csv", append = T)
# 		}
# 	}
# )




# Run simulations (parallel) -------------------------------

run.sim <- TRUE

if (run.sim) {
	library(foreach)
	library(doSNOW)

	# Setup parallel backend to use many processors
	cores <- parallel::detectCores()
	cl <- makeCluster(cores[1]-1)
	registerDoSNOW(cl)

	# Progress bar
	p.bar <- txtProgressBar(max = len.run, style = 3)
	progress <- function(n) setTxtProgressBar(p.bar, n)
	opts <- list(progress = progress)

	alfR::lok_regar( {
		results.df.par <- foreach(i = first.run:max.run, .combine = rbind, .options.snow = opts) %dopar% {
			# Read parameters from table
			gamma <- as.numeric(df.clean[i, "gamma"])
			R1.init <- as.integer(as.numeric(df.clean[i, "R1.init"]) * scale)
			R2.init <- as.integer(as.numeric(df.clean[i, "R2.init"]) * scale)

			# Update parameters
			params <- c(max.time,
									R1.init, R2.init, 0, 0, 0, 0,
									kappa1, omega, gamma, sigma1,
									epsilon1, delta1, curr.capacity,
									n.sim)

			# Run simulation
			prob <- as.numeric(system(paste("./compiled/6eq_prob_exit", paste(params, collapse = ' ')), intern = T))
			results <- cbind(gamma, R1.init, R2.init, prob)
			readr::write_csv(data.frame(results), "data/6eq_ssa_prob_extinction_fucking_big.csv", append = T)
		}
		close(p.bar)
	})

	stopCluster(cl)

	# Reshape results as data frame
	results.df.par <- results.df.par %>%
		data.frame()

	# Save results into file
	# write_csv(results.df.par, "data/6eq_ssa_prob_extinction.csv")
}

# Get results and plot them --------------------------------

# results.df <- data.table::fread("data/6eq_ssa_prob_extinction_fucking_big.csv") %>%
# 	rename(gamma = V1, R1.init = V2, R2.init = V3, prob = V4)
#
# ggplot(results.df %>% filter(R2.init != 0, R1.init != 0, R2.init %in% seq(50,1000,50))) +
# 	aes(x = gamma, y = prob, color = as.factor(R2.init)) +
# 	geom_point() +
# 	geom_line(linetype = "dashed") +
# 	labs(color = "R2.init") +
# 	scale_x_continuous(breaks = seq(0, 0.55, 0.05)) +
# 	scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, 0.1)) +
# 	labs(x = "Gamma", y = "Probability of extinction", color = "Initial R2") +
# 	theme_bw() +
# 	facet_wrap(~R2.init)

