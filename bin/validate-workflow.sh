#!/bin/bash
jq . < "$1" > /dev/null;
ec=$?
if (( ec > 0 )); then
	echo "::error file=$1,line=0,col=0:: Invalid JSON"
	exit 1
fi

topic=$(echo "$1" | cut -d/ -f 2)

jq '.tags[]' -r < "$1" | grep -q "$topic"
ec=$?

if (( ec > 0 )); then
	echo "::warning file=$1,line=0,col=0:: The tags attribute is missing. Please add: \`\"tags\": [\"$topic\"]\`"
	exit_with=1
fi

annotation=$(jq '.annotation' -r < "$1")
if [[ "$annotation" == "null" ]] || [[ "$annotation" == "" ]]; then
	echo "::warning file=$1,line=0,col=0:: The 'annotation' attribute is missing or incorrect. Please add: \`\"annotation\": \"<title of tutorial>\"\`"
	exit_with=1
fi

tts_tools=$(jq '.steps[].tool_id' -r < "$1" | grep -c testtoolshed.g2.bx.psu.edu)
if (( tts_tools > 0 )); then
	echo "::warning file=$1,line=0,col=0:: Some steps have tools from the testtoolshed. These are not permitted in GTN tutorials."
	jq '.steps[].tool_id' -r < "$1" | grep -c testtoolshed.g2.bx.psu.edu
	exit_with=1
fi

exit $exit_with
