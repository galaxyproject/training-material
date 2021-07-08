#!/bin/bash
echo "Checking for broken symlinks"
find -xtype l
count=$(find -xtype l | wc -l)
if (( count > 0 )); then
	echo "Found broken symlinks"
	exit 1
fi
