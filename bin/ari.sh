#!/bin/bash
set -eo pipefail

GTN_CACHE="$(pwd)/.jekyll-cache/aws-polly/"
mkdir -p "$GTN_CACHE"

# Setup our Engine
mozilla_port=$(ss -lptn 'sport = :5002' | wc -l)
if [[ -n "$AWS_ACCESS_KEY_ID" ]]; then
	engine="aws"
	echo "Using AWS for speech (this will cost money!)"
elif (( mozilla_port > 1 )); then
	engine="mozilla"
	echo "Using MozillaTTS for speech"
else
	echo "Cannot run, no speech engine setup. "
	echo "If you're testing, try: docker run -it -p 5002:5002 synesthesiam/mozillatts"
	exit 1
fi

# Setup a temporary directory for building our project.
build_dir=$(mktemp -d /tmp/gtn-ari.XXXXXXXXXX)
echo "Building in $build_dir"

# Setup inputs
if [ "$#" -ne 3 ]; then
    echo "Error, expecting 3 parameters: slides_PDF slides_source mp4_output"
	exit 1
fi

slides=$1  # e.g. _site/training-material/topic/admin/tutorials/ansible/slides.pdf
source=$2  # e.g. topic/admin/tutorials/ansible/slides.html
output=$3  # e.g. _site/training-material/topic/admin/tutorials/ansible/slides.mp4
slidesbase=$(basename "$source" .html)
if [[ "$slidesbase" == *"_ES" ]]; then
	lang=es
	language=esp
else
	lang=en
	language=eng
fi

subtitles="$(dirname "$output")"/"$(basename "$output" .mp4)".${lang}.vtt
cover="$output".png
srcdir="$(dirname "$source")"

# Metadata
meta_authors="$(ruby bin/extract-frontmatter.rb "${source}" | jq '.contributors | join(", ")' -r)"
meta_title="$(ruby bin/extract-frontmatter.rb "${source}" | jq .title -r)"
REVISION="$(git log -1 --format=%H)"

# This is digusting, but the safest way to get a fresh ffmpeg.
FFMPEG_PATH=$(echo "const ffmpeg = require('ffmpeg-static');console.log(ffmpeg.split('/').slice(0, -1).join('/'));" | node -)
echo "Located FFMPEG at $FFMPEG_PATH"
export PATH="$FFMPEG_PATH:$PATH"
which ffmpeg

# We'll cache audio locally.
ffmpeglog=warning

# Setup output files
script="${build_dir}/script.json"

# Generate our script
echo "  Building Script, Subtitles"
echo ruby bin/ari-extract-script.rb "$source"
ruby bin/ari-extract-script.rb "$source" > "$script"
voice=$(cat "$script" | jq .voice.id)

# Now explode that into individual lines, synthesize, and re-assemble them into
# our editly.json5 script
# This will read the environment variable EDITLY_FAST=true and set fast in the script, if requested.
ruby bin/ari-prep-script.rb "${build_dir}" "${engine}"

# Generate images for use.
echo "  Extracting slides"
convert "${slides}" -resize 1920x1080 "${build_dir}/slides.%03d.png"

echo "  Building Video | $(npm bin)/editly --json ${build_dir}/editly.json5"
$FFMPEG_PATH/ffmpeg -version
$(npm bin)/editly --json "${build_dir}/editly.json5"

# Mux it together
echo "  Muxing"
ffmpeg -loglevel $ffmpeglog -i "${build_dir}/tmp.mp4" -i "${build_dir}/out.srt" \
	-movflags +faststart \
	-metadata comment="build-tag:$(date --rfc-3339=seconds)/$REVISION/$USER/$engine/$voice" \
	-metadata network="Galaxy Training Network"\
	-metadata artist="$meta_authors, AWS Polly $voice" \
	-metadata title="$meta_title" \
	-c:v copy -c:a copy -c:s mov_text \
	-map 0 -map 1 \
	-metadata:s:a:0 language=$language \
	-metadata:s:v:0 language=$language \
	-metadata:s:s:0 language=$language \
	-disposition:s:s:0 forced "${build_dir}/out.mp4"

# Check if output dir needs to be created
mkdir -p "$(dirname "$output")"

# Copy our files over
cp "${build_dir}/out.mp4" "$output"
cp "${build_dir}/out.vtt" "$subtitles"
convert -resize 480x270  "${build_dir}/slides.000.png" "$cover"
ffprobe -loglevel warning -show_format -show_private_data -show_streams -print_format json -i "$output" > "$output".json
