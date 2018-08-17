library(data.table)

args <- commandArgs()
# print(args)

file <- args[6]

df <- data.table::fread(paste0(file,"-paint.csv"))
df <- as.data.table(df)

step <-  1/(length(unique(df$V3)) - 1)
count <-  1

surface <- data.frame(matrix(ncol = 3 ))

for (p in seq(min(df$V3), round(max(df$V3), 2), step ) ) {
	if (nrow(df[round(V3, 2) == p, ]) >= 1) {
		for (R2 in seq(min(df$V2), round(max(df$V2), 2), step ) ) {
			R1.vect <- df[round(V2, 2) == R2, ][round(V3, 2) == p, round(V1, 2)]

			if (length(R1.vect) >= 1) {
				surface[count,] <- c(round(max(R1.vect), 2), R2, p)
				count <- count + 1
			}
		}
	}
}

write.table(surface, paste0(file,"-surface.csv"), row.names = F, col.names = F, sep = ",")



