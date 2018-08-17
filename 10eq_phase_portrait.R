# Source: http://tbb.bio.uu.nl/rdb/grindR/

# install.packages(c("deSolve", "rootSolve", "FME"))

source("libs/grind.R")

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

p <- c(kappa1 = 1, kappa2 = 1, kappa3 = 1, kappa4 = 0.5,
			 alpha = 1, omega = 1, mu = 1, beta = 0.6,
			 gamma1 = 0.1, gamma2 = 0.1, gamma3 = 0.1, gamma4 = 0.1,
			 sigma1 = 0.05, sigma2 = 0.05, sigma3 = 0.05,
			 epsilon1 = 0.1, epsilon2 = 0.1, epsilon3 = 0.1,
			 delta1 = 0.05, delta2 = 0.05, delta3 = 0.05) # p is a named vector of parameters


s <- c(RNA1 = 1,
			 RNA2 = 0.1,
			 RNA3 = 0.1,
			 RNA4 = 0,
			 replicase = 0,
			 coat = 0,
			 movement = 0,
			 virion1 = 0,
			 virion2 = 0,
			 virion3 = 0)                # s is the state

# Run simulation
run(tstep = 0.1, tmax = 1000)

# Phase portrait
# pdf(file="phase_pot_fp3_R1_R2_010.pdf",width=6,height=5)

# plane(x = "RNA2", y = "replicase", tstep = 0.5, portrait = T)

# dev.off()

# library(plot3D)
# source("cube.R")
# cube(x="R1", y="R2", z="p1")
#
# run3d(x=1, y=2, z=3, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, tstep=0.1, tmax=10000)
#




