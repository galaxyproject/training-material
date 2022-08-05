#!/bin/bash
ec=0
echo "Checking for broken symlinks"
find -xtype l
count=$(find -xtype l | wc -l)
if (( count > 0 )); then
	echo "Found broken symlinks"
	ec=1
fi

count=$(grep '{%\s*include .*faq' -R topics -c)
if (( count > 0 )); then
	echo "Found snippet issue: you must load FAQs with 'snippet', not 'include'"
	ec=1
fi

exit $ec
