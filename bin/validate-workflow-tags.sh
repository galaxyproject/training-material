#!/bin/bash
exit_with=0

#python script to iterate over the steps in the workflow.
function tester { 
    python3 - <<END
import sys 
import json 
with open("$1") as json_file: 
    data = json.load(json_file) 
    if 'tags' not in data or "$2" not in data['tags']: 
        sys.exit(False)
    else:
        #Checking if there are tools used from the testtoolshed
        for step in data['steps'].values():
            if step['tool_id'] and step['type'] == 'tool' and 'testtoolshed.g2.bx.psu.edu' in step['tool_id']:
                    sys.exit(False)
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
                    echo "ERROR: No 'tags' attribute found with its corresponding topic or workflow includes tool(s) from testtoolshed."
		            exit_with=1
                fi
            done
        fi
    done
done

exit $exit_with
