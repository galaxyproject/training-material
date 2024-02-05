#!/bin/bash
tutorial=topics/dev/tutorials/tool-from-scratch/tutorial.md
cast=$(dirname $tutorial)/$(basename $tutorial .md).cast

# If the tutorial.md is *older* than the tutorial.cast, exit
if [ "$tutorial" -ot "$cast" ]; then
	echo "$tutorial is older than $cast, no need to re-record"
	exit 1
else
	echo "$tutorial is newer than $cast, re-recording"
fi

# Re-record
asciinema rec /tmp/out.cast -t "AUTO_GTN: Tools From Scratch" --command ./bin/video-term-recorder-simple.sh --overwrite --cols 131 --rows 33

# Add chapter markers
python bin/video-postprocess.py /tmp/out.cast topics/dev/tutorials/tool-from-scratch/tutorial.cast
