#!/bin/bash
exit_with=0

function tester { 
    python3 - <<END
import sys 
import json 
with open("$1") as json_file: 
    data = json.load(json_file) 
    if 'tags' not in data or "$2" not in data['tags']: 
        sys.exit(False)
    else:
        sys.exit(True)
END
}

for topicdir in ./topics/*
do
    topic=$(basename $topicdir)
    for tutdir in $topicdir/tutorials/*
    do
        if [ -d $tutdir/workflows/ ];
        then
            for w in $tutdir/workflows/*.ga
            do
                echo "Checking tutorial $w for tags" 
                if tester $w $topic;
                then
                    echo "No 'tags' attribute found in workflow $w with its corresponding topic."
		            exit_with=1
                fi
            done
        fi
    done
done

exit $exit_with
