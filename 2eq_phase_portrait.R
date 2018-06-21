# Source: http://tbb.bio.uu.nl/rdb/grindR/

# install.packages(c("deSolve", "rootSolve", "FME"))

source("libs/grind.R")

model <- function(t, state, parms) {
	with(as.list(c(state,parms)), {
		dR1 <- kappa * R1 * p1 * (1 - R1) - gamma * R1
		dp1 <-  alpha * R1 * (1 - p1) - sigma * p1
		return(list(c(dR1, dp1)))
	})
}

p <- c(kappa = 1,
			 alpha = 1,
			 gamma = 0.1,
			 sigma = 0.1) # p is a named vector of parameters


s <- c(R1 = 1,
			 p1 = 0) # s is the state

# Run simulation
run(tstep = 0.1, tmax = 1000)

# Phase portrait
# pdf(file="phase_pot_fp3_R1_R2_010.pdf",width=6,height=5)

plane(x = "R1", y = "p1", tstep = 0.5, portrait = T)


# dev.off()

# library(plot3D)
# source("cube.R")
# cube(x="R1", y="R2", z="p1")
#
# run3d(x=1, y=2, z=3, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, tstep=0.1, tmax=10000)
#
