directory=data
outfile_name=$3
maxtime=$1
numsim=$2

p0=0.0

# Clean output file
rm $directory/$outfile_name.csv;

# Looping over N
for gamma in 0.00 0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.50 0.51 0.52 0.53 0.54; do
	for (( i = 0; i < 10; i++ )); do

		# Set initial R for each gamma
		if [[ $gamma == 0.00 ]]; then
			R0=1
		elif [[ $gamma == 0.01 ]]; then
			R0=2
		elif [[ $gamma == 0.02 ]]; then
			R0=3
		elif [[ $gamma == 0.03 ]]; then
			R0=4
		elif [[ $gamma == 0.04 ]]; then
			R0=6
		elif [[ $gamma == 0.05 ]]; then
			R0=8
		elif [[ $gamma == 0.10 ]]; then
			R0=20
		elif [[ $gamma == 0.15 ]]; then
			R0=37
		elif [[ $gamma == 0.20 ]]; then
			R0=59
		elif [[ $gamma == 0.25 ]]; then
			R0=87
		elif [[ $gamma == 0.30 ]]; then
			R0=122
		elif [[ $gamma == 0.35 ]]; then
			R0=166
		elif [[ $gamma == 0.40 ]]; then
			R0=222
		elif [[ $gamma == 0.45 ]]; then
			R0=296
		elif [[ $gamma == 0.50 ]]; then
			R0=408
		elif [[ $gamma == 0.51 ]]; then
			R0=444
		elif [[ $gamma == 0.52 ]]; then
			R0=480
		elif [[ $gamma == 0.53 ]]; then
			R0=536
		elif [[ $gamma == 0.54 ]]; then
			R0=1000
		fi

		echo "Gamma: $gamma, iteration: $i, R0: $R0, p0: $p0"

		./compiled/2eq_prob_exit $maxtime $R0 $p0 1 1 $gamma 0.1 $numsim >> $directory/$outfile_name.csv;
	done
done


# 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.50 0.51 0.52 0.53 0.54
