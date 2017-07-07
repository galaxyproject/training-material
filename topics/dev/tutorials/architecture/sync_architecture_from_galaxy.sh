#!/bin/bash

## This script pulls in updates from Galaxy's copy of the the architecture slides and updates this training-material copy.
## Usage: sync_architecture_from_galaxy.sh <path/to/galaxy/root>

GALAXY_ROOT=$1
GALAXY_SLIDES_DIR="$GALAXY_ROOT/doc/source/slideshow/architecture/"
SLIDES_DIR="."

cat > "$SLIDES_DIR/architecture.html" <<EOF
---
layout: tutorial_slides
topic_name: "dev"
tutorial_name: architecture
logo: "GTN"
---
EOF

cat "$GALAXY_ROOT/doc/source/slideshow/architecture/galaxy_architecture.md" >> "$SLIDES_DIR/architecture.html"

awk 'match($0, /images\/.*(svg|png)/) {print substr($0, RSTART, RLENGTH)}' < $GALAXY_SLIDES_DIR/galaxy_architecture.md | xargs -I '{}' cp $GALAXY_SLIDES_DIR/'{}' $SLIDES_DIR/../images
sed 's+(images/+(../../images/+g' "$SLIDES_DIR/architecture.html"
