#!/bin/bash
set -e

# We have an old commit ID, so we need to figure out which slides to build.
if [[ "${PREVIOUS_COMMIT_ID}" != "none" ]]; then
	changed_slides="$(join <(find topics -name 'slides.html' -or -name introduction.html | xargs ./bin/filter-resource-metadata video | sort) <(git diff ${PREVIOUS_COMMIT_ID} --name-only | sort))"
else
	changed_slides="$(find topics -name 'slides.html' -or -name introduction.html | xargs ./bin/filter-resource-metadata video)"
fi

for slides in $changed_slides; do
	echo "====== $slides ======"
	dir="$(dirname "$slides")"
	pdf="$dir/$(basename "$slides" .html).pdf"
	mp4="videos/$dir/$(basename "$slides" .html).mp4"
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
	echo ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
	./bin/ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
done

# Now we'll note our current, changed commit since this all went so well.
mkdir -p videos/topics/
git log -1 --format=%H > videos/topics/last-commit
