#!/bin/bash
for slides in $(find topics -name 'slides.html' | xargs ./bin/filter-has-videos); do
	echo "====== $slides ======"
	dir="$(dirname "$slides")"
	pdf="$dir/$(basename "$slides" .html).pdf"
	mp4="$dir/$(basename "$slides" .html).mp4"
	built_slides="_site/training-material/$slides"

	# Process the slides
	echo $built_slides
	cat "$built_slides" | \
		sed "s|/training-material/|$(pwd)/_site/training-material/|g"  | \
		sed "s|<head>|<head><base href=\"file://$(pwd)/$slides\">|" | \
		wkhtmltopdf \
			--enable-javascript --javascript-delay 3000 --page-width 700px --page-height 530px -B 5px -L 5px -R 5px -T 5px \
			--user-style-sheet bin/slides-fix.css \
			- _site/training-material/"$pdf";

	# Build the slides
	./bin/ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
done
