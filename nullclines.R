nullclines.R <- function(kappa, gamma, num) {
	p = 0;
	R = numeric(num+1);
	p.step = 1/num;
	for(i in 0:num+1){
		R[i] = (kappa * p - gamma) / (kappa * p)
		p = p + p.step;
	}
	return(R)
}

nullclines.p <- function(omega, sigma, num) {
	R = 0;
	p = numeric(num+1);
	R.step = 1/num;
	for(i in 0:num+1) {
		p[i] <- (omega * R) / (omega * R + sigma)
		R = R + R.step;
	}
	return(p)
}

plot_trajectory <- function(kappa, omega, gamma, sigma, step = 1000, scale = 1000) {
	time <- seq(0, 1, 1/step) * scale
	R <- nullclines.R(kappa, gamma, step) * scale
	p <- nullclines.p(omega, sigma, step) * scale
	df <- data.frame(time, R, p)

	gg <- ggplot(df) +
		geom_line(data = df %>% filter(R >= 0), aes(x = R , y = time, colour = "RNA")) +
		geom_line(aes(x = time , y = p, colour = "viral.replicase")) +
		geom_path(data = data.2eq.ssa.401.3, aes(x = V2, y = V3, colour = "path.sto.410")) +
		geom_path(data = data.2eq.ssa.401.4, aes(x = V2, y = V3, colour = "path.sto")) +
		# geom_path(data = data.2eq.ssa.50.800, aes(x = V2, y = V3, colour = "path.sto")) +
		geom_path(data = data.2eq.ssa.400, aes(x = V2, y = V3, colour = "path.sto")) +
		geom_path(data = data.2eq * 1000, aes(x = V2, y = V3, colour = "path.det.410"), linetype = "dashed") +
		geom_path(data = data.2eq.05 * 1000, aes(x = V2, y = V3, colour = "path.det"), linetype = "dashed") +
		# geom_path(data = data.2eq.005.08 * 1000, aes(x = V2, y = V3, colour = "path.det"), linetype = "dashed") +
		geom_path(data = data.2eq.04 * 1000, aes(x = V2, y = V3, colour = "path.det"), linetype = "dashed") +
		labs(x = "RNA", y = "p", color = "Variable") +
		scale_color_manual(values = c("RNA" = "red3",
																	"viral.replicase" = "springgreen4",
																	"path.sto" = "slateblue1",
																	"path.sto.410" = "orchid3",
																	"path.det.410" = "orchid3",
																	"path.det" = "slateblue1")) +
		theme_bw()

	return(gg)
}

gg.trajs <- plot_trajectory(1, 1, 0.5, 0.1)

gg.trajs

saveRDS(gg.trajs, "trajectories.rds")

gg.traj <- readRDS("trajectories.rds")
