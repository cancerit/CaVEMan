########## LICENSE ##########
# Copyright (c) 2014-2015 Genome Research Ltd. 
# 
# Author: Cancer Genome Project cgpit@sanger.ac.uk 
# 
# This file is part of CaVEMan. 
# 
# CaVEMan is free software: you can redistribute it and/or modify it under 
# the terms of the GNU Affero General Public License as published by the Free 
# Software Foundation; either version 3 of the License, or (at your option) any 
# later version. 
# 
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more 
# details. 
# 
# You should have received a copy of the GNU Affero General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
###########################

echo "Running unit tests:"

for i in tests/*_tests
do
	if test -f $i
	then
		if $VALGRIND ./$i 2>> tests/tests_log
		then
			echo $i PASS
		else
			echo "ERROR in test $i: here's tests/tests_log"
			echo "------"
			tail tests/tests_log
			exit 1
		fi
	fi
done

echo ""