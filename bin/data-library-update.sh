#!/bin/bash

for group in $(find topics -name tutorial.md); do
	echo $group
	topic=$(echo "$group" | cut -d/ -f 2)
	tutos=$(echo "$group" | cut -d/ -f 4)

	zenodo=$(ruby bin/extract-frontmatter.rb $group | jq '.zenodo_link | length')
	if (( zenodo > 0 )); then
		planemo training_fill_data_library --topic_name $topic --tutorial_name $tutos
	fi
done
