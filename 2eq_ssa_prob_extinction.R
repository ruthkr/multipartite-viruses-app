library(tidyverse)

source("2eq_find_first_point_to_P2.R")

# Compile
# system("gcc -Wall -Ofast -lm 2eq_ssa_prob_exit.c -o compiled/2eq_prob_exit")

# General parameters
max.time <- 1000
n.sim <- 1000

scale <- 1000
n.run <- nrow(df.clean)

# Run simulations (standard) ----------------------------------------------

# Create data frame
# results <- data.frame(cbind(0,0,0,0))
# results <- results %>%
# 	dplyr::rename("gamma" = X1, "R1.init" = X2, "p.init" = X3, "prob" = X4)
#
# alfR::lok_regar(
# 	for (i in 1:n.run) {
# 		print(paste0("Simulation ", i, "/", n.run, " started at ", format(Sys.time(), "%H:%M:%S")))
# 		# Read parameters from table
# 		gamma <- as.numeric(df.clean[i, "gamma"])
# 		R1.init <- as.integer(as.numeric(df.clean[i, "R1.init"]) * scale)
# 		p.init <- as.integer(as.numeric(df.clean[i, "p.init"]) * scale)
#
# 		# Update parameters
# 		params <- c(max.time,
# 								R1.init, p.init,
# 								kappa, omega, gamma, sigma,
# 								n.sim)
#
# 		# Run simulation
# 		prob <- as.numeric(system(paste("./compiled/2eq_prob_exit", paste(params, collapse = ' ')), intern = T))
# 		results[1,] <- cbind(gamma, R1.init, p.init, prob)
#
# 		# print(results)
# 		if (i == 1) {
# 			write_csv(results, "data/2eq_ssa_prob_extinction_old.csv")
# 		} else {
# 			write_csv(results, "data/2eq_ssa_prob_extinction_old.csv", append = T)
# 		}
# 	}
# )


# Run simulations (parallel) -------------------------------
run.sim <- FALSE

if (run.sim) {
	library(foreach)
	library(doSNOW)

	# Setup parallel backend to use many processors
	cores <- parallel::detectCores()
	cl <- makeCluster(cores[1]-1)
	registerDoSNOW(cl)

	# Progress bar
	p.bar <- txtProgressBar(max = n.run, style = 3)
	progress <- function(n) setTxtProgressBar(p.bar, n)
	opts <- list(progress = progress)

	alfR::lok_regar({
		results.df.par <- foreach(i = 1:n.run, .combine = rbind, .options.snow = opts) %dopar% {
			# Read parameters from table
			gamma <- as.numeric(df.clean[i, "gamma"])
			R1.init <- as.integer(as.numeric(df.clean[i, "R1.init"]) * scale)
			p.init <- as.integer(as.numeric(df.clean[i, "p.init"]) * scale)

			# Update parameters
			params <- c(max.time,
									R1.init, p.init,
									kappa, omega, gamma, sigma,
									n.sim)

			# Run simulation
			prob <- as.numeric(system(paste("./compiled/2eq_prob_exit", paste(params, collapse = ' ')), intern = T))
			cbind(gamma, R1.init, p.init, prob)
		}
		close(p.bar)
	})

	stopCluster(cl)

	# Reshape results as data frame
	results.df.par <- results.df.par %>%
		data.frame()

	# Save results into file
	# write_csv(results.df.par, "data/2eq_ssa_prob_extinction.csv")
}

# Get results and plot them --------------------------------

results.df <- data.table::fread("data/2eq_ssa_prob_extinction.csv")

gg.2eq.extinction <- ggplot(results.df) +
	aes(x = gamma, y = prob, color = as.factor(p.init)) +
	geom_point() +
	geom_line(linetype = "dashed") +
	scale_x_continuous(breaks = seq(0, 1, 0.1)) +
	labs( x = "Gamma", y = "Probabilty of extinction", color = "Initial replicase") +
	theme_bw()

gg.2eq.extinction + theme(text = element_text(family = "LM Roman 10")) + ggsave(filename = "2eq_prob_extinction.pdf", width = 6, height = 3.5, dpi = 96, device = "pdf")

gg.2eq.extinction.small <- ggplot(results.df) +
	aes(x = gamma, y = prob, color = as.factor(p.init)) +
	geom_point(size = 0.5) +
	geom_line(linetype = "dashed", size = 0.5) +
	scale_x_continuous(breaks = seq(0, 1, 0.2)) +
	labs( x = "Gamma", y = "Probabilty of extinction", color = "Initial replicase") +
	theme_bw()

gg.2eq.extinction.facet <- gg.2eq.extinction.small + facet_wrap(~ p.init, ncol = 3)

gg.2eq.extinction.facet + theme(strip.background = element_blank(), strip.text.x = element_blank()) + theme(text = element_text(family = "LM Roman 10")) + ggsave(filename = "2eq_prob_extinction_facet.pdf", width = 6, height = 6, dpi = 96, device = "pdf")

