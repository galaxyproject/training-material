#!/bin/bash
USED=$(mktemp)
ALL=$(mktemp)


for topic in $(find topics -maxdepth 1 -mindepth 1 -type d); do
	echo $topic

	# get directly referenced images
	egrep '!\[[^[]*\]\([^)]*\)' $topic/tutorials/*/slides.html $topic/tutorials/*/tutorial.md -h 2>/dev/null| \
		grep ']([^(]*)' -o | \
		grep -v '://' | \
		sed 's/^](//' | \
		grep -F .. | \
		sed 's/ .*//g;s/)$//' | \
		sed "s|\.\./\.\.|$topic|" | \
		sed "s|)$||g" | \
		xargs -n1 -I {} realpath --relative-base=$(pwd) {} \
		>> $USED.all

	# Linked
	egrep '!\[[^[]*\]\([^)]*\)' $topic/tutorials/*/slides.html $topic/tutorials/*/tutorial.md -h 2>/dev/null | \
		grep ']([^(]*)' -o | \
		grep -v '://' | \
		sed 's/^](//' | \
		grep -v -F .. | \
		grep topics | \
		sed 's/{{\s*site.baseurl\s*}}{% link //g;s/\s*%}.*//g;s/{{site.baseurl}}\///g;s/ .*//' | \
		sed "s|)$||g" \
		>> $USED.all


	# Intro slides are special
	# get directly referenced images
	egrep '!\[[^[]*\]\([^)]*\)' $topic/slides/*.html -h 2>/dev/null | \
		grep ']([^(]*)' -o | \
		grep -v '://' | \
		sed 's/^](//' | \
		grep -F .. | \
		sed 's/ .*//g;s/)$//' | \
		sed "s|\.\.|$topic|" | \
		sed "s|)$||g" | \
		xargs -n1 -I {} realpath --relative-base=$(pwd) {} \
		>> $USED.all

	# Linked
	egrep '!\[[^[]*\]\([^)]*\)' $topic/slides/*.html -h 2>/dev/null | \
		grep ']([^(]*)' -o | \
		grep -v '://' | \
		sed 's/^](//' | \
		grep -v -F .. | \
		grep topics | \
		sed 's/{{\s*site.baseurl\s*}}{% link //g;s/\s*%}.*//g;s/{{site.baseurl}}\///g;s/ .*//' | \
		sed "s|)$||g" \
		>> $USED.all
done

find topics/*/images/ -type f > $ALL.all
find shared/images -type f >> $ALL.all

sort -u $USED.all > $USED
sort -u $ALL.all > $ALL
rm -f $USED.all $ALL.all

echo "Missing images"
diff $USED $ALL | grep '<'

echo "Unusued images"
diff $USED $ALL | grep '>' | sed 's/> //g'| xargs -I{} rm '{}'

rm -f $USED $ALL
