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
else
	lang=en
fi

subtitles="$(dirname "$output")"/"$(basename "$output" .mp4)".${lang}.vtt
cover="$output".png
srcdir="$(dirname "$source")"

# Metadata
meta_authors="$(ruby bin/extract-frontmatter.rb topics/admin/tutorials/cvmfs/slides.html | jq '.contributors | join(", ")' -r)"
meta_title="$(ruby bin/extract-frontmatter.rb "${source}" | jq .title -r)"
REVISION="$(git log -1 --format=%H)"

# We'll cache audio locally.
ffmpeglog=warning

# Setup output files
script="${build_dir}/script.json"
images="${build_dir}/images.txt"
sounds="${build_dir}/sounds.txt"

# Generate our script
echo ruby bin/ari-extract-script.rb "$source"
ruby bin/ari-extract-script.rb "$source" > "$script"

# Now explode that into individual lines, synthesize, and re-assemble them into
# our images.txt/sounds.txt scripts
ruby bin/ari-prep-script.rb "${script}" "${build_dir}" "${engine}"


# Generate images for use.
echo "  Extracting slides"
convert "${slides}" -resize 1920x1080 "${build_dir}/slides.%03d.png"

# Generate a pause of 1 seconds, used after slides.
sox -n -r 44100 -c 2 "${build_dir}/silence-unfixed.mp3" trim 0 1
ffmpeg -loglevel $ffmpeglog -i "${build_dir}/silence-unfixed.mp3"  "${build_dir}/silence.mp3"

# Build components
echo "  Building Subtitles"
python3 bin/ari-subs.py "${build_dir}" --format srt > "${build_dir}/tmp.srt"
python3 bin/ari-subs.py "${build_dir}" --format webvtt > "${build_dir}/out.vtt"
echo "  Building Audio"


# Concatenate (literally) the files. This is how the concat filter (not
# demuxer!!!) works in ffmpeg. It makes NO sense. but it gets the right output time.
# If you use the concat demuxer, you get something that's off by ~3 seconds
# over 250 seconds. It's not clear what generates this discrepancy but it's
# solved by first concatenating them all, then transcoding once within the same
# format (still wrong time) then transcoding again to the new format (finally
# correct.)
for f in $(cat "${build_dir}/sounds.txt" | grep '^file' | cut -f2 -d' ' | sed "s/'//g"); do
	cat "${build_dir}/$f" >> "${build_dir}/ugly.mp3"
done
# This conversion step fixes the headers (still the wrong timestamps)
ffmpeg -loglevel fatal -i "${build_dir}/ugly.mp3" "${build_dir}/ugly2.mp3"

# Finally it's right.
ffmpeg -loglevel fatal -i "${build_dir}/ugly2.mp3" "${build_dir}/tmp.m4a"
echo "  Building Video"
ffmpeg -loglevel $ffmpeglog -f concat -i "$images" -pix_fmt yuv420p -vcodec h264 "${build_dir}/tmp.mp4"

# Mux it together
echo "  Muxing"
ffmpeg -loglevel $ffmpeglog -i "${build_dir}/tmp.mp4" -i "${build_dir}/tmp.m4a" -i "${build_dir}/tmp.srt" \
	-movflags +faststart \
	-metadata comment="build-tag:$(date --rfc-3339=seconds)/$REVISION/$USER/$engine" \
	-metadata network="Galaxy Training Network"\
	-metadata artist="$meta_authors" \
	-metadata title="$meta_title" \
	-c:v copy -c:a copy -c:s mov_text \
	-map 0:v:0 -map 1:a:0 -map 2 \
	-metadata:s:s:0 language=English -disposition:s:s:0 forced \
	-b:a 192k "${build_dir}/out.mp4"

# Check if output dir needs to be created
mkdir -p "$(dirname "$output")"

# Copy our files over
cp "${build_dir}/out.mp4" "$output"
cp "${build_dir}/out.vtt" "$subtitles"
cp "${build_dir}/slides.000.png" "$cover"
ffprobe -loglevel warning -show_format -show_private_data -show_streams -print_format json -i "$output" > "$output".json
