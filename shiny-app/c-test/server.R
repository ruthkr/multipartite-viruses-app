#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

	# Initial conditions

	R1.init <- 1
	R2.init <- 0.4
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
	n.iter <- 100

	# Flags
	fixed.point.select <- "third"

	params <- c(n.iter,
							R1.init, R2.init, p.init, s.init,
							v1.init, v2.init,
							kappa1, kappa2, alpha, beta, gamma1, gamma2, sigma1, sigma2,
							epsilon1, epsilon2, delta1, delta2)

	#Compile
	system("gcc -Ofast -lm rk4_bipartite_6eq_test.c -o rk4-test")
	system(paste("./rk4-test", paste(params, collapse = ' '), "> results-test.csv"))

	# Read data
	data <- data.table::fread("results-test.csv")

  output$distPlot <- renderPlot({
  	ggplot(data) +
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
  })

})
