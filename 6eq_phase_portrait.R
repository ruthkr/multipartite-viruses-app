# Source: http://tbb.bio.uu.nl/rdb/grindR/

# install.packages(c("deSolve", "rootSolve", "FME"))

source("libs/grind.R")

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

p <- c(kappa1 = 1, kappa2 = 1,
			 alpha = 1, beta = 1,
			 gamma1 = 0.2, gamma2 = 0.2,
			 sigma1 = 0.1, sigma2 = 0.1,
			 epsilon1 = 0.1, epsilon2 = 0.1,
			 delta1 = 0.1, delta2 = 0.1) # p is a named vector of parameters


s <- c(R1 = 1,
			 R2 = 0.07,
			 p1 = 0,
			 s1 = 0,
			 v1 = 0,
			 v2 = 0)                # s is the state

# Run simulation
# run(tstep = 0.1, tmax = 1000)

# Phase portrait
# pdf(file="phase_pot_fp3_R1_R2_010.pdf",width=6,height=5)

plane(x = "R1", y = "R2", tstep = 0.5, portrait = T)


# dev.off()

 # library(plot3D)
 # source("cube.R")
 # cube(x="R1", y="R2", z="p1")
 #
 # run3d(x=1, y=2, z=3, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, tstep=0.1, tmax=10000)
 #




