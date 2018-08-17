# Source: http://tbb.bio.uu.nl/rdb/grindR/

# install.packages(c("deSolve", "rootSolve", "FME"))

source("libs/grind.R")

model <- function(t, state, parms) {
	with(as.list(c(state,parms)), {
		dRNA1 <- kappa1 * RNA1 * replicase * (1 - RNA1- RNA2) - gamma1 * RNA1 - epsilon1 * RNA1 * coat
		dRNA2 <- kappa2 * RNA2 * replicase * (1 - RNA1- RNA2) - gamma2 * RNA2 - epsilon2 * RNA2 * coat
		dreplicase <-  alpha * RNA1 * (1 - replicase - coat) - sigma1 * replicase
		dcoat <-  beta * RNA2 * (1 - replicase - coat) - sigma2 * coat
		dvirion1 <- epsilon1 * RNA1 * coat - delta1 * virion1
		dvirion2 <- epsilon2 * RNA2 * coat - delta2 * virion2
		return(list(c(dRNA1, dRNA2, dreplicase, dcoat, dvirion1, dvirion2)))
	})
}

p <- c(kappa1 = 1, kappa2 = 1,
			 alpha = 1, beta = 1,
			 gamma1 = 0.5, gamma2 = 0.5,
			 sigma1 = 0.1, sigma2 = 0.1,
			 epsilon1 = 0.1, epsilon2 = 0.1,
			 delta1 = 0.1, delta2 = 0.1) # p is a named vector of parameters


s <- c(RNA1 = 1,
			 RNA2 = 0.7,
			 replicase = 0,
			 coat = 0,
			 virion1 = 0,
			 virion2 = 0)                # s is the state

# Run simulation
# run(tstep = 0.1, tmax = 1000)

# Phase portrait
# pdf(file="phase_pot_fp3_R1_R2_010.pdf",width=6,height=5)

plane(x = "RNA1", y = "RNA2", tstep = 0.5, portrait = T)

# dev.off()

 # library(plot3D)
 # source("cube.R")
 # cube(x="R1", y="R2", z="p1")
 #
 # run3d(x=1, y=2, z=3, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, tstep=0.1, tmax=10000)
 #




