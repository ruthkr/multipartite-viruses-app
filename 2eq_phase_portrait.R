# Source: http://tbb.bio.uu.nl/rdb/grindR/

# install.packages(c("deSolve", "rootSolve", "FME"))

source("libs/grind.R")

model <- function(t, state, parms) {
	with(as.list(c(state,parms)), {
		dRNA1 <- kappa * RNA1 * replicase * (1 - RNA1) - gamma * RNA1
		dreplicase <-  alpha * RNA1 * (1 - replicase) - sigma * replicase
		return(list(c(dRNA1, dreplicase)))
	})
}

p <- c(kappa = 0.5,
			 alpha = 1,
			 gamma = 0.1,
			 sigma = 0.1) # p is a named vector of parameters


s <- c(RNA1 = 1,
			 replicase = 0) # s is the state

# Run simulation
# run(tstep = 0.1, tmax = 1000)

# Phase portrait
# pdf(file="phase_pot_fp3_RNA1_R2_010.pdf",width=6,height=5)

plane(x = "RNA1", y = "replicase", tstep = 0.5, portrait = T)


# dev.off()

# library(plot3D)
# source("cube.R")
# cube(x="RNA1", y="R2", z="replicase")
#
# run3d(x=1, y=2, z=3, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, tstep=0.1, tmax=10000)
#
