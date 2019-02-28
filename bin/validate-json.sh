#!/bin/bash
exit_with=0

# For each worfklow
for wf in $(find topics -name '*.ga'); do
	# Check that it is valid JSON
	echo "Checking $wf"
	jq . < $wf > /dev/null;
	ec=$?
	if (( ec > 0 )); then
		echo "  Invalid JSON"
		exit_with=1
	fi
done

exit $exit_with
