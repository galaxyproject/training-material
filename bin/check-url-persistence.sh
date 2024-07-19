#!/bin/bash
if [[ ! -f /tmp/2021.txt ]]; then
	curl --silent https://training.galaxyproject.org/archive/2021-02-01/sitemap.xml | sed 's|<url>|\n|g' | grep '<loc>[^<]*</loc>' -o | sed 's/<loc>//;s/<\/loc>//' | sed 's|archive/2021-02-01|training-material|g'  > /tmp/2021.txt
fi

if [[ ! -f /tmp/2022.txt ]]; then
	curl --silent https://training.galaxyproject.org/archive/2022-01-01/sitemap.xml | sed 's|<url>|\n|g' | grep '<loc>[^<]*</loc>' -o | sed 's/<loc>//;s/<\/loc>//' | sed 's|archive/2022-01-01|training-material|g'  > /tmp/2022.txt
fi

if [[ ! -f /tmp/2024.txt ]]; then
	curl --silent https://training.galaxyproject.org/archive/2024-01-01/sitemap.xml | sed 's|<url>|\n|g' | grep '<loc>[^<]*</loc>' -o | sed 's/<loc>//;s/<\/loc>//' | sed 's|archive/2024-01-01|training-material|g'  > /tmp/2024.txt
fi


if [[ ! -f /tmp/2099.txt ]]; then
	curl --silent https://training.galaxyproject.org/training-material/sitemap.xml | sed 's|<url>|\n|g' | grep '<loc>[^<]*</loc>' -o | sed 's/<loc>//;s/<\/loc>//' > /tmp/2099.txt
fi

# No guarantees of API or data-file persistence
# 1fe4d7d92e5ea5a5794cbe741eadb96a74511261
cat /tmp/20*.txt | sort -u | \
	grep -v '/api/' | grep -v '/by-tool/'  | grep -v '/hall-of-fame/' | \
	grep -v '/badges/' | \
	grep --extended-regexp -v 'krona_?[a-z]*.html' | \
	grep -v '/transcriptomics/tutorials/ref-based/faqs/rnaseq_data.html' | \
	grep -v '/topics/data-management/' | \
	grep -v 'training-material/tags/' | grep -v 'data-library'| \
	sed 's|/$|/index.html|'  | grep '.html$' | sort -u | sed 's|https://training.galaxyproject.org|_site|' > /tmp/gtn-files.txt

count=0
while read line; do
	if [[ ! -f $line ]]; then
		echo "Missing: $line"
		count=$((count+1))
	fi
done < /tmp/gtn-files.txt

if (( $count > 0 )); then
	echo "Files in previous versions that are not currently redirected: $count"
	echo "Please ensure you add redirect_from: entries to each page"
	echo "If this file is intentionally deleted, with no replacement, please update this script with some exclusions."
	exit 1
fi
