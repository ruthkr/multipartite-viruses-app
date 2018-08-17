library(tidyverse)

source("2eq_find_first_point_to_P2.R")

scale <- 1000
gg_2eq_separatrix_point <- ggplot(df.clean %>%
																		filter(gamma %in% c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5) || gamma == 0.53) %>%
																		mutate(R1.init = R1.init * scale,
																					 p.init = p.init * scale)) +
	geom_point(aes(x = (R1.init + 2) , y = p.init, color = as.factor(gamma) ) ) +
	geom_line(aes(x = (R1.init - 1), y = p.init, colour = as.factor(gamma)), linetype = "dashed") +
	labs(x = "RNA1", y = "Replicase", color = "Gamma") +
	scale_color_discrete(labels = c("0.00", "0.05", "0.10", "0.15",
																	"0.20", "0.25", "0.30", "0.35",
																	"0.40", "0.45", "0.50", "0.53")) +
	scale_x_continuous(breaks = seq(0,1000,100)) +
	theme_bw()

gg_2eq_separatrix_point #+ theme(text = element_text(family = "LM Roman 10")) + ggsave(filename = "gg_2eq_separatrix_point.pdf", width = 6, height = 3.5, dpi = 96, device = "pdf")
# saveRDS(gg_2eq_separatrix_point, paste0("objects/", "gg_2eq_separatrix_point.rds"))
