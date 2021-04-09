#!/usr/bin/env bash
set -e

# A script that performs all steps needed to upgrade tours.
if [[ $# != 1 ]] ; then
  echo "USAGE: $0 /output/directory"
  exit 1
fi
set -u

if [[ ! -e $1 ]] ; then
  echo "$1 does not exist. Creating it."
  mkdir -p $1
fi

python bin/fix_tours.py $1

for tour in $(find $1 -name tour.yaml.j2) ; do
  dir=$(dirname $tour)
  echo "Generating tour $tour"
  python bin/update_tour.py $dir
  python bin/render_tour.py $dir
done



