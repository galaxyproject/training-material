#!/bin/bash
output=$(mktemp)

# Check the toolshed-y ones.
grep -n '\{% tool \[' -R topics | sed 's/:\([0-9]\+\):.*{% tool/:\1/g;s/)\s*%}.*//g;s/ \[/\t/g;s/](/\t/' | \
	awk -F'\t' '(substr($3,1,18) != "toolshed.g2.bx.psu"){ print $0 }' | \
	egrep -i -v '(Count1|Cut1|Filter1|Grep1|Grouping1|Remove\+beginning1|Show\+beginning1|Summary_Statistics1|addValue|cat1|comp1|gene2exon1|join1|random_lines1|ucsc_table_direct1|upload1)' | \
	awk -F'\t' '{print $1"\t"$3}' | tee $output

# If there was output, exit.
lines=$(wc -l $output | cut -f1 -d' ')
if (( lines > 0 )); then
	echo "ERROR: There were multiple broken tool links. Please see above and correct."
	exit 1;
fi

rm -f $output
