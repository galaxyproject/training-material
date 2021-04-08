#!/usr/bin/env bash

if [[ $# != 1 ]] ; then
  echo "USAGE: $0 /output/directory"
  exit 1
fi

if [[ ! -e $1 ]] ; then
  echo "Could not find $1"
  exit 1
fi

python bin/fix_tours.py $1

for tour in $(find $1 -name tour.yaml.j2) ; do
  dir=$(dirname $tour)
  echo "Generating tour $tour"
  python bin/update_tour.py $dir
  python bin/render_tour.py $dir
done



