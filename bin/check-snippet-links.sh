#!/bin/bash
IFS=$'\n'
ec=0
for row in $(grep -n '{% snippet [^ ]* .*%}' -R topics | sed -e 's/:\([0-9]*\):.*{% snippet /\t\1\t/;s/ .*//'); do
	fn=$(echo "$row" | cut -f 3)
	if [[ ! -f "$fn" ]]; then
		echo "$row" | awk '{print "::error file="$1",line="$2",col=0:: The referenced snippet file does not exist." }'
		ec=1
	fi
done
exit $ec

