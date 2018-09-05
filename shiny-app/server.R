# Libraries ------------------------------------------------

# Server Libs
# library(dplyr)
# library(xts)
library(shinyjs)

# Server ---------------------------------------------------

server <- function(input, output) {
	source("libs/grind.R") # For phase portrait
	source("libs/cube.R")  # For 3D trajectories

	# Dynamic plot list
	# observeEvent(input$run_sims, {
	output$ui_plots <- renderUI({
		out <- list()
		if (length(input$`plots_checkbox`)==0){return(NULL)}

		for (i in input$`plots_checkbox`){
			# The name of the repl_* renderPlot()s must match the checkboxGroupInput()s
			out[[i]] <- plotOutput(outputId = paste0(input$sidebarmenu, "_", i))
		}
		return(out)
		# })
	})

	# TS Replication ####
	output$replication_time_series <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa * RNA1 * replicase * (1 - RNA1) - gamma * RNA1
				dp1 <-  alpha * RNA1 * (1 - replicase) - sigma * replicase
				return(list(c(dR1, dp1)))
			})
		}

		p <- c(kappa = input$slider_kappa,
					 alpha = input$slider_alpha,
					 gamma = input$slider_gamma,
					 sigma = input$slider_sigma)


		s <- c(RNA1 = input$slider_r1,
					 replicase = 0)

		run(tstep = 0.5, state = s, odes = model, parms = p)
	})

	# PP Replication ####
	output$replication_phase_port <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa * RNA1 * replicase * (1 - RNA1) - gamma * RNA1
				dp1 <-  alpha * RNA1 * (1 - replicase) - sigma * replicase
				return(list(c(dR1, dp1)))
			})
		}

		p <- c(kappa = input$slider_kappa,
					 alpha = input$slider_alpha,
					 gamma = input$slider_gamma,
					 sigma = input$slider_sigma)


		s <- c(RNA1 = input$slider_r1,
					 replicase = 0)

		plane(x = "RNA1", y = "replicase", tstep = 0.5, state = s, odes = model, parms = p, portrait = T)
	})

	# TS Bipartite ####
	output$bipartite_time_series <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
				dR2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
				dp1 <-  alpha * RNA1 * (1 - replicase - coat) - sigma1 * replicase
				ds1 <-  beta * RNA2 * (1 - replicase - coat) - sigma2 * coat
				dv1 <- epsilon1 * RNA1 * coat - delta1 * virion1
				dv2 <- epsilon2 * RNA2 * coat - delta2 * virion2
				return(list(c(dR1, dR2, dp1, ds1, dv1, dv2)))
			})
		}

		p <- c(
			kappa1   = input$slider_kappa_bi,   # default: 1
			kappa2   = input$slider_kappa_bi,
			alpha    = input$slider_alpha_bi,   # default: 1
			beta     = input$slider_alpha_bi,
			gamma1   = input$slider_gamma_bi,   # default: 0.1
			gamma2   = input$slider_gamma_bi,
			sigma1   = input$slider_sigma_bi,   # default: 0.1
			sigma2   = input$slider_sigma_bi,
			epsilon1 = input$slider_epsilon_bi, # default: 0.1
			epsilon2 = input$slider_epsilon_bi,
			delta1   = input$slider_delta_bi,   # default: 0.1
			delta2   = input$slider_delta_bi
		)

		s <- c(
			RNA1 = input$slider_r1_bi,
			RNA2 = input$slider_r2_bi,
			replicase = 0,
			coat = 0,
			virion1 = 0,
			virion2 = 0
		)

		run(tstep = 0.5, state = s, odes = model, parms = p)
	})

	# Traj Bipartite ####
	output$bipartite_3d_traj <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
				dR2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
				dp1 <-  alpha * RNA1 * (1 - replicase - coat) - sigma1 * replicase
				ds1 <-  beta * RNA2 * (1 - replicase - coat) - sigma2 * coat
				dv1 <- epsilon1 * RNA1 * coat - delta1 * virion1
				dv2 <- epsilon2 * RNA2 * coat - delta2 * virion2
				return(list(c(dR1, dR2, dp1, ds1, dv1, dv2)))
			})
		}

		p <- c(
			kappa1   = input$slider_kappa_bi,   # default: 1
			kappa2   = input$slider_kappa_bi,
			alpha    = input$slider_alpha_bi,   # default: 1
			beta     = input$slider_alpha_bi,
			gamma1   = input$slider_gamma_bi,   # default: 0.1
			gamma2   = input$slider_gamma_bi,
			sigma1   = input$slider_sigma_bi,   # default: 0.1
			sigma2   = input$slider_sigma_bi,
			epsilon1 = input$slider_epsilon_bi, # default: 0.1
			epsilon2 = input$slider_epsilon_bi,
			delta1   = input$slider_delta_bi,   # default: 0.1
			delta2   = input$slider_delta_bi
		)

		s <- c(
			RNA1 = input$slider_r1_bi,
			RNA2 = input$slider_r2_bi,
			replicase = 0,
			coat = 0,
			virion1 = 0,
			virion2 = 0
		)

		cube(x = 1, y = 2, z = 3, state = s)
	})

	# PP Bipartite ####
	output$bipartite_phase_port <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
				dR2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
				dp1 <-  alpha * RNA1 * (1 - replicase - coat) - sigma1 * replicase
				ds1 <-  beta * RNA2 * (1 - replicase - coat) - sigma2 * coat
				dv1 <- epsilon1 * RNA1 * coat - delta1 * virion1
				dv2 <- epsilon2 * RNA2 * coat - delta2 * virion2
				return(list(c(dR1, dR2, dp1, ds1, dv1, dv2)))
			})
		}

		p <- c(
			kappa1   = input$slider_kappa_bi,   # default: 1
			kappa2   = input$slider_kappa_bi,
			alpha    = input$slider_alpha_bi,   # default: 1
			beta     = input$slider_alpha_bi,
			gamma1   = input$slider_gamma_bi,   # default: 0.1
			gamma2   = input$slider_gamma_bi,
			sigma1   = input$slider_sigma_bi,   # default: 0.1
			sigma2   = input$slider_sigma_bi,
			epsilon1 = input$slider_epsilon_bi, # default: 0.1
			epsilon2 = input$slider_epsilon_bi,
			delta1   = input$slider_delta_bi,   # default: 0.1
			delta2   = input$slider_delta_bi
		)

		s <- c(
			RNA1 = input$slider_r1_bi,
			RNA2 = input$slider_r2_bi,
			replicase = 0,
			coat = 0,
			virion1 = 0,
			virion2 = 0
		)

		plane(x = "RNA1", y = "RNA2", tstep = 0.5, state = s, odes = model, parms = p, portrait = T)
	})

	# TS Tripartite ####
	output$tripartite_time_series <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dRNA1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
				dRNA2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
				dRNA3 <- kappa3 * RNA3 * replicase * (1 - beta) * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma3 * RNA3 - epsilon3 * RNA3 * coat
				dRNA4 <- kappa4 * RNA3 * replicase * beta * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma4 * RNA4
				dreplicase <-  alpha^2 * RNA1 * RNA2 * (1 - replicase - coat - movement) - sigma1 * replicase
				dcoat <-  omega * RNA4 * (1 - replicase - coat - movement) - sigma2 * coat
				dmovement <- mu * RNA3 * (1 - replicase - coat - movement) - sigma3 * movement
				dvirion1 <- epsilon1 * RNA1 * coat - delta1 * virion1 - movement * virion1
				dvirion2 <- epsilon2 * RNA2 * coat - delta2 * virion2 - movement * virion2
				dvirion3 <- epsilon3 * RNA3 * coat - delta3 * virion3 - movement * virion3
				return(list(c(dRNA1, dRNA2, dRNA3, dRNA4, dreplicase, dcoat, dmovement, dvirion1, dvirion2, dvirion3)))
			})
		}

		p <- c(
			kappa1 = input$slider_kappa_tri,
			kappa2 = input$slider_kappa_tri,
			kappa3 = input$slider_kappa_tri,
			kappa4 = input$slider_kappa_tri,
			alpha = input$slider_alpha_tri,
			omega = input$slider_omega_tri,
			mu = input$slider_mu_tri,
			beta = input$slider_beta_tri,
			gamma1 = input$slider_gamma_tri,
			gamma2 = input$slider_gamma_tri,
			gamma3 = input$slider_gamma_tri,
			gamma4 = input$slider_gamma_tri,
			sigma1 = input$slider_sigma_tri,
			sigma2 = input$slider_sigma_tri,
			sigma3 = input$slider_sigma_tri,
			epsilon1 = input$slider_epsilon_tri,
			epsilon2 = input$slider_epsilon_tri,
			epsilon3 = input$slider_epsilon_tri,
			delta1 = input$slider_delta_tri,
			delta2 = input$slider_delta_tri,
			delta3 = input$slider_delta_tri
		)


		s <- c(
			RNA1 = input$slider_r1_tri,
			RNA2 = input$slider_r2_tri,
			RNA3 = input$slider_r3_tri,
			RNA4 = 0,
			replicase = 0,
			coat = 0,
			movement = 0,
			virion1 = 0,
			virion2 = 0,
			virion3 = 0
		)

		# Run simulation
		run(tstep = 0.1, tmax = 1000, state = s, odes = model, parms = p)
	})

	# PP Tripartite  ####
	output$tripartite_phase_port <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dRNA1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
				dRNA2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
				dRNA3 <- kappa3 * RNA3 * replicase * (1 - beta) * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma3 * RNA3 - epsilon3 * RNA3 * coat
				dRNA4 <- kappa4 * RNA3 * replicase * beta * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma4 * RNA4
				dreplicase <-  alpha^2 * RNA1 * RNA2 * (1 - replicase - coat - movement) - sigma1 * replicase
				dcoat <-  omega * RNA4 * (1 - replicase - coat - movement) - sigma2 * coat
				dmovement <- mu * RNA3 * (1 - replicase - coat - movement) - sigma3 * movement
				dvirion1 <- epsilon1 * RNA1 * coat - delta1 * virion1 - movement * virion1
				dvirion2 <- epsilon2 * RNA2 * coat - delta2 * virion2 - movement * virion2
				dvirion3 <- epsilon3 * RNA3 * coat - delta3 * virion3 - movement * virion3
				return(list(c(dRNA1, dRNA2, dRNA3, dRNA4, dreplicase, dcoat, dmovement, dvirion1, dvirion2, dvirion3)))
			})
		}

		p <- c(
			kappa1 = input$slider_kappa_tri,
			kappa2 = input$slider_kappa_tri,
			kappa3 = input$slider_kappa_tri,
			kappa4 = input$slider_kappa_tri,
			alpha = input$slider_alpha_tri,
			omega = input$slider_omega_tri,
			mu = input$slider_mu_tri,
			beta = input$slider_beta_tri,
			gamma1 = input$slider_gamma_tri,
			gamma2 = input$slider_gamma_tri,
			gamma3 = input$slider_gamma_tri,
			gamma4 = input$slider_gamma_tri,
			sigma1 = input$slider_sigma_tri,
			sigma2 = input$slider_sigma_tri,
			sigma3 = input$slider_sigma_tri,
			epsilon1 = input$slider_epsilon_tri,
			epsilon2 = input$slider_epsilon_tri,
			epsilon3 = input$slider_epsilon_tri,
			delta1 = input$slider_delta_tri,
			delta2 = input$slider_delta_tri,
			delta3 = input$slider_delta_tri
		)


		s <- c(
			RNA1 = input$slider_r1_tri,
			RNA2 = input$slider_r2_tri,
			RNA3 = input$slider_r3_tri,
			RNA4 = 0,
			replicase = 0,
			coat = 0,
			movement = 0,
			virion1 = 0,
			virion2 = 0,
			virion3 = 0
		)

		plane(x = "RNA1", y = "RNA2", tstep = 0.1, state = s, odes = model, parms = p, portrait = T)
	})

	# Traj Tripartite  ####
	output$tripartite_3d_traj <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dRNA1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
				dRNA2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
				dRNA3 <- kappa3 * RNA3 * replicase * (1 - beta) * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma3 * RNA3 - epsilon3 * RNA3 * coat
				dRNA4 <- kappa4 * RNA3 * replicase * beta * (1 - RNA1- RNA2 - RNA3 - RNA4) - gamma4 * RNA4
				dreplicase <-  alpha^2 * RNA1 * RNA2 * (1 - replicase - coat - movement) - sigma1 * replicase
				dcoat <-  omega * RNA4 * (1 - replicase - coat - movement) - sigma2 * coat
				dmovement <- mu * RNA3 * (1 - replicase - coat - movement) - sigma3 * movement
				dvirion1 <- epsilon1 * RNA1 * coat - delta1 * virion1 - movement * virion1
				dvirion2 <- epsilon2 * RNA2 * coat - delta2 * virion2 - movement * virion2
				dvirion3 <- epsilon3 * RNA3 * coat - delta3 * virion3 - movement * virion3
				return(list(c(dRNA1, dRNA2, dRNA3, dRNA4, dreplicase, dcoat, dmovement, dvirion1, dvirion2, dvirion3)))
			})
		}

		p <- c(
			kappa1 = input$slider_kappa_tri,
			kappa2 = input$slider_kappa_tri,
			kappa3 = input$slider_kappa_tri,
			kappa4 = input$slider_kappa_tri,
			alpha = input$slider_alpha_tri,
			omega = input$slider_omega_tri,
			mu = input$slider_mu_tri,
			beta = input$slider_beta_tri,
			gamma1 = input$slider_gamma_tri,
			gamma2 = input$slider_gamma_tri,
			gamma3 = input$slider_gamma_tri,
			gamma4 = input$slider_gamma_tri,
			sigma1 = input$slider_sigma_tri,
			sigma2 = input$slider_sigma_tri,
			sigma3 = input$slider_sigma_tri,
			epsilon1 = input$slider_epsilon_tri,
			epsilon2 = input$slider_epsilon_tri,
			epsilon3 = input$slider_epsilon_tri,
			delta1 = input$slider_delta_tri,
			delta2 = input$slider_delta_tri,
			delta3 = input$slider_delta_tri
		)


		s <- c(
			RNA1 = input$slider_r1_tri,
			RNA2 = input$slider_r2_tri,
			RNA3 = input$slider_r3_tri,
			RNA4 = 0,
			replicase = 0,
			coat = 0,
			movement = 0,
			virion1 = 0,
			virion2 = 0,
			virion3 = 0
		)

		cube(x = 1, y = 2, z = 3, state = s)
	})


	# Debugger
	output$debugger <- renderText({
		paste(input$plots_checkbox)
	})

}
