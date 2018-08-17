#!/usr/bin/awk -f

# Authors: Ruth Kristianingsih <ruth.kristianingsih@mathmods.eu>

# Legal Stuff:
#	This script is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, version 3.

#	This script is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#	GNU General Public License for more details.

#	You should have received a copy of the GNU General Public License
#	along with this script. If not, see <http://www.gnu.org/licenses/>.

BEGIN{
	FS = ", "
	OFS = ","
}
{
	if ($1 == gamma) {
		print $5, $6, $7
	}
}

# INSTRUCTIONS:
# awk -v gamma=0.0 -f clean_3d_data.awk mathematica_gamma_3d.csv

# && ($4 == 0 || $4 == 0.75 || $4 == 0.6 || $4 == 0.9
