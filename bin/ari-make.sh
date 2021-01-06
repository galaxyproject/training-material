#!/bin/bash
for slides in $(find topics -name 'slides.html' | xargs ./bin/filter-has-videos); do
	echo "====== $slides ======"
	dir="$(dirname "$slides")"
	pdf="$dir/$(basename "$slides" .html).pdf"
	mp4="$dir/$(basename "$slides" .html).mp4"

	# Build the slides
	make "_site/training-material/$pdf" ACTIVATE_ENV=pwd

	# Build the slides
	./bin/ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
done
