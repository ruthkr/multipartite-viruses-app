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
	select(gamma, count, ratio)


g <- ggplot(data=a, aes(x=gamma, y=ratio, group=1))  +
	geom_line(linetype = "dashed") +
	labs(x = "Gamma", y = "Ratio") +
	geom_point(size = 2, colour = "springgreen4") +
	theme_bw()

g + theme(axis.text.x=element_text(angle=50, vjust=0.5))
