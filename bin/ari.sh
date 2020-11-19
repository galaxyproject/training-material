#!/bin/bash
VOICE_ID=Amy
POLLY_PARAMS="--engine neural --language-code en-GB --voice-id $VOICE_ID"
GTN_CACHE="$(pwd)/.jekyll-cache/aws-polly/"
mkdir -p "$GTN_CACHE"

# Setup a temporary directory for building our project.
build_dir=$(mktemp -d /tmp/gtn-ari.XXXXXXXXXX)
echo "Building in $build_dir"

# Setup inputs
slides=$1; shift; # _site/training-material/topic/admin/tutorials/ansible/slides.pdf
source=$1; shift; # topic/admin/tutorials/ansible/slides.html
output=$1; shift; # _site/training-material/topic/admin/tutorials/ansible/slides.mp4
aws=$1; shift; # Empty or 'upload'
srcdir="$(dirname "$source")"

# Extract the script from a tutorial
generate_script(){
	# The title page
	./bin/extract-frontmatter.rb "$1" | jq .title -r
	# If there's a questions page, generate a line reading "Questions"
	question_count=$(./bin/extract-frontmatter.rb "$1" | jq '.questions | length')
	if (( question_count > 0 )); then
		./bin/extract-frontmatter.rb "$1" | jq '.questions[]' -r | paste -s -d' '
	fi

	# If there's a objectives page, generate a line reading "Objectives"
	objective_count=$(./bin/extract-frontmatter.rb "$1" | jq '.objectives | length')
	if (( objective_count > 0 )); then
		./bin/extract-frontmatter.rb "$1" | jq '.objectives[]' -r | paste -s -d' '
	fi

	# The contents
	cat "$1" | \
		sed -n '/???/p;/???/,/---/{//!p};' | \
		tr '\n' ' ' | \
		sed 's/^??? //' | \
		sed 's/???\s*/\n/g' | \
		sed 's/\s*$/   /g' | \
		pandoc -t html | \
		sed -r 's|</?p>||g;s|<br />$||g' | \
		sed 's/<[^>]*>/ /g'
}

cache_speech() {
	line="$1"; shift;

	linehash="$(echo "$line" | md5sum | cut -c1-32 | sed 's/\s*//g')"

	fn="${GTN_CACHE}/${linehash}-${VOICE_ID}"
	# Generate audio
	if [[ ! -f "$fn.mp3" ]]; then
		aws polly synthesize-speech $POLLY_PARAMS --output-format mp3 --text "$line" "$fn.mp3" &>/dev/null
	fi

	# And sentence markers
	if [[ ! -f "$fn.json" ]]; then
		aws polly synthesize-speech $POLLY_PARAMS --output-format json --text "$line" --speech-mark-types='["sentence"]' "$fn.json" &>/dev/null
	fi

	# Copy from the cache to our build directory
	cp "$fn.mp3"  "${build_dir}/${linehash}.mp3"
	cp "$fn.json" "${build_dir}/${linehash}.json"

	echo $linehash
}

# We'll cache audio locally.
ffmpeglog=warning

# Setup output files
script="${build_dir}/script.txt"
images="${build_dir}/images.txt"
sounds="${build_dir}/sounds.txt"

# Generate our script
generate_script "$source" > "$script"

echo "  There are $(wc -w "$script" | awk '{print $1}') words in the script"

i=0;
while read -r line; do
	# Process this line of the script

	echo "  Speaking/$VOICE_ID: $(echo "$line" | cut -c1-64)..."
	# Create the speech file
	hash=$(cache_speech "$line" "${build_dir}" | sed 's/\s*//g')

	# Get the duration of the clip
	duration=$(ffprobe -show_streams -print_format json "${build_dir}/${hash}.mp3" 2>/dev/null | jq .streams[0].duration -r)

	# Write out our slides. Original slide + extra second pause
	printf "file 'slides.%03d.png'\nduration %s\n" $i "$duration" >> "$images"
	printf "file 'slides.%03d.png'\nduration 1.04\n" $i >> "$images"

	# Write out our audio, original + 1 second of silence.
	printf "file '%s.mp3'\nduration %s\n" "$hash" "$duration" >> "$sounds"
	printf "file 'silence.mp3'\nduration 1.04\n" >> "$sounds"

	i=$((i + 1));
done < "$script"

# Generate images for use.
echo "  Extracting slides"
convert -density 300 "${slides}" "${build_dir}/slides.%03d.png"

# Thank you for the end slides
hash=$(cache_speech 'Thank you for watching!' "${build_dir}")
duration=$(ffprobe -show_streams -print_format json "${build_dir}/${hash}.mp3" 2>/dev/null | jq .streams[0].duration -r)
last_slide=$(basename $(find $build_dir -name '*.png' | sort | tail -n 1))

printf "file '%s'\nduration %s\n" $last_slide "$duration" >> "$images"
printf "file '%s'\n" $last_slide >> "$images" # Bug in concat filter.
printf "file '%s.mp3'\nduration %s\n" "$hash" "$duration" >> "$sounds"

# Generate a pause of 1 seconds, used after slides.
sox -n -r 44100 -c 2 ${build_dir}/silence.mp3 trim 0 1

# Build components
echo "  Building Subtitles"
python3 bin/ari-subs.py "${build_dir}" --format srt > "${build_dir}/tmp.srt"
python3 bin/ari-subs.py "${build_dir}" --format webvtt > "${build_dir}/out.vtt"
echo "  Building Audio"
# Concat first within the same format, this is very important. If you don't, the time codes get massively screwed up
ffmpeg -loglevel $ffmpeglog -f concat -i "$sounds" "${build_dir}/tmp.mp3"
# Then convert to the format we really want, m4a.
ffmpeg -loglevel $ffmpeglog -i "${build_dir}/tmp.mp3" "${build_dir}/tmp.m4a"
echo "  Building Video"
ffmpeg -loglevel $ffmpeglog -f concat -i "$images" -pix_fmt yuv420p -vcodec h264 "${build_dir}/tmp.mp4"
# Mux it together
echo "  Muxing"
ffmpeg -loglevel $ffmpeglog -i "${build_dir}/tmp.mp4" -i "${build_dir}/tmp.m4a" -i "${build_dir}/tmp.srt" \
	-c:v copy -c:a copy -c:s mov_text -map 0:v:0 -map 1:a:0 -map 2 -b:a 192k "${build_dir}/out.mp4"

cp "${build_dir}/out.mp4" "$output"

if [[ "$aws" == "upload" ]]; then
	aws s3 cp "${build_dir}/out.mp4" "s3://galaxy-training/videos/$srcdir/slides.mp4"
	aws s3 cp "${build_dir}/out.vtt" "s3://galaxy-training/videos/$srcdir/slides.en.vtt"
fi
