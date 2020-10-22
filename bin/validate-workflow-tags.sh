#!/bin/bash
exit_with=0

for wf in $(find topics -path '*/workflows/*' -and -name '*.ga'); do
	topic=$(echo "$wf" | cut -d/ -f 2)

	echo "Linting $wf"

	jq '.tags[]' -r < "$wf" | grep -q "$topic"
	ec=$?

	if (( ec > 0 )); then
		echo "	The tags attribute is missing. Please add:"
		echo '	"tags": ["'"$topic"'"]'
		exit_with=1
	fi

	annotation=$(jq '.annotation' -r < "$wf")
	if [[ "$annotation" == "null" ]] || [[ "$annotation" == "" ]]; then
		echo "	The 'annotation' attribute is missing or incorrect. Please add: \"annotation\": \"<title of tutorial>\""
		exit_with=1
	fi

	tts_tools=$(jq '.steps[].tool_id' -r < "$wf" | grep -c testtoolshed.g2.bx.psu.edu)
	if (( tts_tools > 0 )); then
		echo "	Some steps have tools from the testtoolshed. These are not permitted in GTN tutorials."
		jq '.steps[].tool_id' -r < "$wf" | grep -c testtoolshed.g2.bx.psu.edu
		exit_with=1
	fi

done

exit $exit_with