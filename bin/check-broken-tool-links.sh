#!/bin/bash
output=$(mktemp)

# Check the toolshed-y ones.
grep -n '\{% tool \[' -R topics | sed 's/:\([0-9]\+\):.*{% tool/\t\1/g;s/)\s*%}.*//g;s/ \[/\t/g;s/](/\t/' | \
	awk -F'\t' '(substr($4,1,18) != "toolshed.g2.bx.psu"){ print $0 }' | \
	egrep -i -v '(intermine|Count1|Cut1|Filter1|Grep1|Grouping1|Paste1|Remove\+beginning1|Show\+beginning1|Summary_Statistics1|addValue|cat1|comp1|gene2exon1|join1|random_lines1|ucsc_table_direct1|upload1|sort1|wig_to_bigWig|ChangeCase)' | \
    egrep -v '__[A-Z_]+__' | \
	egrep -v 'interactive_tool_[a-z]+' | \
	awk -F'\t' '{print "::warning file="$1",line="$2",col=0::Broken tool link, please double check."}' | tee $output

# If there was output, exit.
lines=$(wc -l $output | cut -f1 -d' ')
if (( lines > 0 )); then
	(>&2 echo "ERROR: There were multiple broken tool links. Please see above and correct.")
	exit 1;
fi

rm -f $output
