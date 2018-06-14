# Libraries ------------------------------------------------

# Server Libs
# library(dplyr)
# library(xts)
library(shinyjs)

# Server ---------------------------------------------------

server <- function(input, output) {
	source("libs/grind.R") # For phase portrait

	# Dynamic plot list
	observeEvent(input$run_sims, {
		output$ui_plots <- renderUI({
			out <- list()
			if (length(input$`plots_checkbox`)==0){return(NULL)}

			for (i in input$`plots_checkbox`){
				# The name of the repl_* renderPlot()s must match the checkboxGroupInput()s
				out[[i]] <- plotOutput(outputId = paste0(input$sidebarmenu, "_", i))
			}
			return(out)
		})
	})

	# Replication model
	output$replication_time_series <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa * R1 * p1 * (1 - R1) - gamma * R1
				dp1 <-  alpha * R1 * (1 - p1) - sigma * p1
				return(list(c(dR1, dp1)))
			})
		}

		p <- c(kappa = input$slider_kappa,
					 alpha = input$slider_alpha,
					 gamma = input$slider_gamma1,
					 sigma = input$slider_sigma)


		s <- c(R1 = input$slider_r1,
					 p1 = 0)

		run(x = "R1", y = "p1", tstep = 0.5, state = s, odes = model, parms = p)
	})

	output$replication_phase_port <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa * R1 * p1 * (1 - R1) - gamma * R1
				dp1 <-  alpha * R1 * (1 - p1) - sigma * p1
				return(list(c(dR1, dp1)))
			})
		}

		p <- c(kappa = input$slider_kappa,
					 alpha = input$slider_alpha,
					 gamma = input$slider_gamma1,
					 sigma = input$slider_sigma)


		s <- c(R1 = input$slider_r1,
					 p1 = 0)

		plane(x = "R1", y = "p1", tstep = 0.5, state = s, odes = model, parms = p, portrait = T)
	})

	# Bipartite model
	output$bipartite_time_series <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa1 * R1 * p1 * (1 - R1- R2) - gamma1 * R1 - epsilon1 * R1 * s1
				dR2 <- kappa2 * R2 * p1 * (1 - R1- R2) - gamma2 * R2 - epsilon2 * R2 * s1
				dp1 <-  alpha * R1 * (1 - p1 - s1) - sigma1 * p1
				ds1 <-  beta * R2 * (1 - p1 - s1) - sigma2 * s1
				dv1 <- epsilon1 * R1 * s1 - delta1 * v1
				dv2 <- epsilon2 * R2 * s1 - delta2 * v2
				return(list(c(dR1, dR2, dp1, ds1, dv1, dv2)))
			})
		}

		p <- c(kappa1 = input$slider_kappa, kappa2 = 1,
					 alpha = 1, beta = 1,
					 gamma1 = input$slider_gamma1, gamma2 = input$slider_gamma1,
					 sigma1 = 0.6, sigma2 = 0.1,
					 epsilon1 = 0.1, epsilon2 = 0.1,
					 delta1 = 0.1, delta2 = 0.1)

		s <- c(R1 = input$slider_r1,
					 R2 = 0.4,
					 p1 = 0.5,
					 s1 = 0,
					 v1 = 0,
					 v2 = 0)

		run(x = "R1", y = "R2", tstep = 0.5, state = s, odes = model, parms = p, portrait = T)
	})

	output$bipartite_phase_port <- renderPlot({
		model <- function(t, state, parms) {
			with(as.list(c(state,parms)), {
				dR1 <- kappa1 * R1 * p1 * (1 - R1- R2) - gamma1 * R1 - epsilon1 * R1 * s1
				dR2 <- kappa2 * R2 * p1 * (1 - R1- R2) - gamma2 * R2 - epsilon2 * R2 * s1
				dp1 <-  alpha * R1 * (1 - p1 - s1) - sigma1 * p1
				ds1 <-  beta * R2 * (1 - p1 - s1) - sigma2 * s1
				dv1 <- epsilon1 * R1 * s1 - delta1 * v1
				dv2 <- epsilon2 * R2 * s1 - delta2 * v2
				return(list(c(dR1, dR2, dp1, ds1, dv1, dv2)))
			})
		}

		p <- c(kappa1 = input$slider_kappa, kappa2 = 1,
					 alpha = 1, beta = 1,
					 gamma1 = input$slider_gamma1, gamma2 = input$slider_gamma1,
					 sigma1 = 0.6, sigma2 = 0.1,
					 epsilon1 = 0.1, epsilon2 = 0.1,
					 delta1 = 0.1, delta2 = 0.1)

		s <- c(R1 = input$slider_r1,
					 R2 = 0.4,
					 p1 = 0.5,
					 s1 = 0,
					 v1 = 0,
					 v2 = 0)

		plane(x = "R1", y = "R2", tstep = 0.5, state = s, odes = model, parms = p, portrait = T)
	})


	# Debugger
	output$debugger <- renderText({
		paste(input$plots_checkbox)
	})

}
