generate_script(){
	./bin/extract-frontmatter.rb $1 | jq .title -r
	cat $1 | \
		sed -n '/???/p;/???/,/---/{//!p};' | \
		tr '\n' ' ' | \
		sed 's/^??? //' | \
		sed 's/???\s*/\n/g' | \
		pandoc -t html | \
		sed -r 's|</?p>||g;s|<br />$||g'
}

# Figure out the Inputs
pdf="_site/training-material/topics/admin/tutorials/ansible/slides.pdf"
html="_site/training-material/topics/admin/tutorials/ansible/slides.html"
dir="_site/training-material/topics/admin/tutorials/ansible/"
source="topics/admin/tutorials/ansible/slides.html"

# Now we need to get the script
script=$(mktemp)
generate_script $source > $script

images=images.txt
sounds=sounds.txt

i=0;

while read -r line; do
	# Process this line of the script

	# Hash the lines so we don't download unnecessarily
	linehash="$(echo "$line" | md5sum | cut -c1-32)"
	# Create it if it doesn't exist
	if [[ ! -f "$linehash.mp3" ]]; then
		aws polly synthesize-speech --engine standard --language-code en-US --output-format mp3 --text "$line" --voice-id Gwyneth $linehash.mp3
	fi

	# Get the duration of the clip
	duration=$(ffprobe -show_streams -print_format json $linehash.mp3 2>/dev/null | jq .streams[0].duration -r)
	# Add a tiny bit
	duration=$(echo "$duration + 0.5" | bc -l)

	# Write out the relevant information to our file
	printf "file 'slides.%03d.png'\nduration %s\n" $i $duration >> $images
	printf "file '%s.mp3'\noutpoint %s\n" $linehash $duration   >> $sounds

	i=$((i + 1));
done < $script

echo ffmpeg -f concat -safe 0 -i $images -f concat -safe 0 -i $sounds -r 25 -pix_fmt yuv420p out.mp4
ffmpeg -f concat -safe 0 -i $images -f concat -safe 0 -i $sounds -r 25 -pix_fmt yuv420p out.mp4
