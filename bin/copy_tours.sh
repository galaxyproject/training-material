#!/usr/bin/env bash
set -eu

# Copy tour files to a local Galaxy installation. Tours will be renamed using
# the name of parent directory.

if [[ $# != 2 ]] ; then
  echo "USAGE: $0 /tours/directory /galaxy/tours/directory"
  exit 1
fi

for tour in $(find $1 -name tour.yaml) ; do
  dir=$(dirname $tour)
  tour_name=$(basename $dir)
  echo "cp $tour $2/$tour_name.yaml"
  cp $tour $2/$tour_name.yaml
done