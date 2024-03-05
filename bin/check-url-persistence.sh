#!/bin/bash
if [[ ! -f /tmp/2022.txt ]]; then
	curl --silent https://training.galaxyproject.org/archive/2022-01-01/sitemap.xml | sed 's|<url>|\n|g' | grep '<loc>[^<]*</loc>' -o | sed 's/<loc>//;s/<\/loc>//' | sed 's|archive/2022-01-01|training-material|g'  > /tmp/2022.txt
fi

if [[ ! -f /tmp/2024.txt ]]; then
	curl --silent https://training.galaxyproject.org/archive/2024-01-01/sitemap.xml | sed 's|<url>|\n|g' | grep '<loc>[^<]*</loc>' -o | sed 's/<loc>//;s/<\/loc>//' | sed 's|archive/2024-01-01|training-material|g'  > /tmp/2024.txt
fi

# No guarantees of API or data-file persistence
cat /tmp/202*.txt | grep -v '/api/' | grep -v '/by-tool/'  | grep -v '/hall-of-fame/.*/index.html' | grep -v 'training-material/tags/' | grep -v 'data-library'| sed 's|/$|/index.html|'  | grep '.html$' | sort -u | sed 's|https://training.galaxyproject.org|_site|' > /tmp/gtn-files.txt

count=0
while read line; do
	if [[ ! -f $line ]]; then
		echo "Missing: $line"
		count=$((count+1))
	fi
done < /tmp/gtn-files.txt

if (( $count > 0 )); then
	echo "Files in previous versions that are not currently redirected: $count"
	exit 1
fi
