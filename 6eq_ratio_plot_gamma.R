library(tidyverse)

a <- system("wc -l data/sim100*-paint.csv", intern = T) %>%
	data.frame() %>%
	head(-1) %>%
	rename("raw" = ".") %>%
	mutate(
		count = as.numeric(substring(raw, 1, 8)),
		gamma = as.factor(substring(raw, 21, 24)),
		gamma = sub('^\\.|\\-$', '', gamma),
		ratio = count/(101^3)
	) %>%
	select(gamma, count, ratio) %>%
	mutate(gamma = as.numeric(gamma))

a.mod <- rbind(
	a,
	c(0.54, 100, 1.0)
)

g <- ggplot(data=a.mod, aes(x=gamma, y=ratio, group=1))  +
	geom_line(linetype = "dashed") +
	geom_point(data = a, size = 2, colour = "springgreen4") +
	geom_vline(xintercept = 0.54, linetype = "dashed", color = "dodgerblue2") +
	scale_x_continuous(breaks = seq(0, 0.7, 0.05)) +
	labs(x = "Gamma", y = expression(paste("Fraction of Initial Conditions Achieving ", P[1]^"*")) ) +
	theme(axis.text.x=element_text(angle=50, vjust=0.5)) +
	theme_bw()


g

saveRDS(g, "objects/plot6eq_ratio_gamma.rds")





