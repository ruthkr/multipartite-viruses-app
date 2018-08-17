library(tidyverse)

source("6eq_find_first_point_to_P3.R")

scale <- 1000

df.gg.sep <- df.clean.sep %>%
	filter(gamma %in% c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53)) %>%
	mutate(R1.init = R1.init * scale,
				 R2.init = R2.init * scale) %>%
	group_by(gamma, R1.init) %>%
	summarise(R2.init = min(R2.init))


df.gg.init <- df.clean.init %>%
	filter(gamma %in% c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.53)) %>%
	mutate(R1.init = R1.init * scale,
				 R2.init = R2.init * scale)

gg_6eq_separatrix_point <- ggplot(df.gg.init) +
	geom_point(aes(x = (R1.init + 2) , y = R2.init, color = as.factor(gamma) ) ) +
	geom_line(data = df.gg.sep, aes(x = (R1.init - 1), y = R2.init, colour = as.factor(gamma)), linetype = "dashed") +
	labs(x = "RNA1", y = "RNA2", color = "Gamma") +
	scale_color_discrete(labels = c("0.00", "0.01", "0.02", "0.03", "0.04", "0.05", "0.10", "0.15",
																	"0.20", "0.25", "0.30", "0.35",
																	"0.40", "0.45", "0.50", "0.53")) +
	scale_x_continuous(breaks = seq(0,1000,100)) +
	scale_y_continuous(breaks = seq(0,1000,100)) +
	theme_bw()

gg_6eq_separatrix_point #+ theme(text = element_text(family = "LM Roman 10")) + ggsave(filename = "gg_2eq_separatrix_point.pdf", width = 6, height = 3.5, dpi = 96, device = "pdf")
# saveRDS(gg_6eq_separatrix_point, paste0("objects/", "gg_2eq_separatrix_point_parR2_10.rds"))
