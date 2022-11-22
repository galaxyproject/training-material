#!/bin/bash
set -e

function cleanup(){
	kill $(pgrep -f $(npm bin)/http-server) || true
}

trap cleanup EXIT

# hack to just do ours for the moment
# eg set BUILD_ME_SLIDES using env to "topics/galaxy-project/slides/introduction.html"
$(npm bin)/http-server -p 9876 _site &
 
for slides in $BUILD_ME_SLIDES; do
	echo "====== $slides ======"
	dir="$(dirname "$slides")"
	pdf="$dir/$(basename "$slides" .html).pdf"
	mp4="videos/$dir/$(basename "$slides" .html).mp4"
	built_slides="_site/training-material/$slides"

	# Process the slides
	echo $built_slides
    docker run --rm --network host -v $(pwd):/slides astefanutti/decktape automatic -s 1920x1080 http://127.0.0.1:9876/training-material/$slides /slides/_site/training-material/$pdf

	# Build the slides
	echo ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
	./bin/ari.sh "_site/training-material/$pdf" "$slides" "$mp4"
done

cleanup


# Now we'll note our current, changed commit since this all went so well.
mkdir -p videos/topics/
git log -1 --format=%H > videos/topics/last-commit
