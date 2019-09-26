#!/bin/bash
ALL=$(mktemp)
echo "Locating images ($ALL)"
find topics/*/images/ -type f > $ALL

overallec=0

for img in $(cat "$ALL"); do
	image_name=$(basename "$img")
	grep -RI "$image_name" topics > /dev/null
	ec=$?

	if (( ec != 0 )); then
		if [[ $1 == "--remove" ]]; then
			echo "Removing $img"
			rm "$img"
		else
			echo "$img appears to be unused"
			overallec=1
		fi
	fi
done

exit $overallec
