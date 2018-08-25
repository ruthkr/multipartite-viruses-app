RDStoCSV <- function(basename, path = "data/", cols = NULL, header = T) {
	RDS.file <- paste0(path, basename, ".rds")
	CSV.file <- paste0(path, gsub("\\.", "_", basename), ".csv")

	RDS.df <- readRDS(RDS.file)

	if (is.numeric(cols)) {
		RDS.df <- RDS.df %>%
			select(cols)
	}

	readr::write_csv(RDS.df, CSV.file, col_names = header)
}


#Example how to run it
RDStoCSV("6eq.det.R1080.R21.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.det.R1064.R2080.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.det.R1048.R2060.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.det.R1032.R2040.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.det.R1016.R2020.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.sto.R1160.R2200.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.sto.R1320.R2400.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.sto.R1480.R2600.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.sto.R1640.R2800.gamma02", cols = c(2,3,4), header = F)
RDStoCSV("6eq.sto.R1800.R21000.gamma02", cols = c(2,3,4), header = F)
