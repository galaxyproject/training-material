#!/bin/bash
set -e

function cleanup(){
	kill $(pgrep -f $(npm bin)/http-server) || true
}

trap cleanup EXIT

# We have an old commit ID, so we need to figure out which slides to build.
videos="$(find topics -name 'slides.html' -or -name introduction.html -or -name 'slides_*ES.html')"
if [[ "${PREVIOUS_COMMIT_ID}" != "none" ]]; then
	changed_slides="$(join <(echo "$videos" | xargs ./bin/filter-resource-metadata video | sort) <(git diff ${PREVIOUS_COMMIT_ID} --name-only | sort))"
else
	changed_slides="$(echo "$videos" | xargs ./bin/filter-resource-metadata video)"
fi

$(npm bin)/http-server -p 9876 _site &

for slides in $changed_slides; do
	echo "====== $slides ======"
	dir="$(dirname "$slides")"
	pdf="$dir/$(basename "$slides" .html).pdf"
	mp4="videos/$dir/$(basename "$slides" .html).mp4"
	built_slides="_site/training-material/$slides"

	# Process the slides
	echo $built_slides
	docker run --rm --network host -v $(pwd):/slides astefanutti/decktape  automatic -s 1920x1080 http://127.0.0.1:9876/training-material/$slides /slides/_site/training-material/$pdf

	# Build the slides
	echo ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
	./bin/ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
done

cleanup


# Now we'll note our current, changed commit since this all went so well.
mkdir -p videos/topics/
git log -1 --format=%H > videos/topics/last-commit
