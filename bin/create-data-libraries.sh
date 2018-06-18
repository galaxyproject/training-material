#!/usr/bin/env bash

for tuto_dir in topics/*/tutorials/*
do
    tuto="$(basename $tuto_dir)"
    topic="$(basename $(dirname $(dirname $tuto_dir)))"
    echo "$topic / $tuto"
    python bin/create-data-library.py --tutorial $tuto --topic $topic
    echo ""
done
