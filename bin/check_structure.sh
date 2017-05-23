#!/usr/bin/env bash
#
# Check if each tutorial have the required files

set -e

function check_dir {
    parent_dir=$1
    dir=$2
    if [ ! -d "$parent_dir/$dir" ]; then
        echo "$dir directory is missing for $parent_dir"
        #exit 1
    fi
}

function check_file {
    parent_dir=$1
    filename=$2
    if [[ ! -f "$parent_dir/$filename" || ! -s "$parent_dir/$filename" ]]; then
        echo "The $filename file is missing or empty for $parent_dir"
        #exit 1
    fi
}

for topic in $(find topics/* -maxdepth 0 -type d); do
    echo $topic
    # Check metadata
    check_file $topic "metadata.yaml"
    # Check README
    check_file $topic "README.md"
    # Check images
    check_dir $topic "images"
    # Check Docker
    check_dir $topic "docker"
    check_file $topic "docker/Dockerfile"
    # Check slides
    check_dir $topic "slides"
    check_file $topic "slides/index.html"
    # Check tutorials
    check_dir $topic "tutorials"
    for tutorial in $(find $topic/tutorials/* -maxdepth 0 -type d); do
        echo $tutorial
        check_file "$tutorial" "tutorial.md"
        check_file "$tutorial" "metadata.yaml"
        check_file "$tutorial" "tools.yaml"
        check_file "$tutorial" "data-library.yaml"
        check_file "$tutorial" "workflow.ga"
        check_file "$tutorial" "tour.yaml"
    done
    echo ""
done