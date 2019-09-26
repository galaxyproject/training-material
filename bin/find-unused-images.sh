#!/bin/bash
echo "Locating images"

overallec=0

while read -r img; do
	image_name="$(basename "$img")"
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
done < <( find topics/*/images/ -type f )

exit $overallec
