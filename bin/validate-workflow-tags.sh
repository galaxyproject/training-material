#!/bin/bash
exit_with=0

#python script to iterate over the steps in the workflow.
function tester { 
    python3 - <<END
import sys 
import json 

problems = 0
output = ["-----------------------ERROR-----------------------------"]
with open("$1") as json_file:
    data = json.load(json_file)
    # Checking for 'tags' in workflow if topic is known
    if 'tags' not in data or not data['tags'] or "$2" not in data['tags']:
        problems += 1
        output.append(
            "{}. The 'tags' attribute is missing. Please add:".format(str(problems), data['name']))
        output.append('"tags": [' + "\n\t" + '"' + "$2" + '"' + "\n]")

    # Checking for 'annotation' in workflow
    if 'annotation' not in data or not data['annotation']:
        problems += 1
        output.append(
            "{}. The 'annotation' attribute is missing. Please add:".format(str(problems)))
        output.append('"annotation": "<title of tutorial>"')

    # Checking if there are tools used from the testtoolshed
    for stepnr, step in data['steps'].items():
        if step['tool_id'] and step['type'] == 'tool' and 'testtoolshed.g2.bx.psu.edu' in step['tool_id']:
            problems += 1
            output.append("{}. Step {} has a tool from the testtoolshed.".format(str(problems), str(stepnr)))

    if problems:
        output.insert(1, "Workflow '{}' has {} problem(s) because:".format(data['name'], str(problems)))
        output.append("---------------------------------------------------------\n")
        sys.stderr.write("\n".join(output))
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
                echo "Checking $w"

                if tester $w $topic;
                then
		            exit_with=1
                fi
            done
        fi
    done
done

exit $exit_with