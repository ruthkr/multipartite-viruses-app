RDStoCSV <- function(basename, path = "data/") {
	RDS.file <- paste0(path, basename, ".rds")
	CSV.file <- paste0(path, gsub("\\.", "_", basename), ".csv")
	readr::write_csv(readRDS(RDS.file), CSV.file)
}


#Example how to run it
RDStoCSV("6eq.sto.R1800.R21000.gamma02")
