#!/bin/bash
SCRIPT="$1";
SCRIPT_LOG="$2";
# Disable git's paging which can cause a hang
export GIT_PAGER=cat
# Prevent pip from shouting everywhere
pip config --user set global.progress_bar off
# Setup the demo-magic
. ./demo-magic/demo-magic.sh -n
# Go into right dir
cd /home/hxr/arbeit/galaxy/git-gat/
# Hide our demo-magic activation
clear

OLDIFS="$IFS"
IFS=$'\n' # bash specific
scriptStart=$(date +%s.%N)
for line in $(cat $SCRIPT | jq '.[]' -c); do
	# These aren't fast.
	action=$(echo "$line" | jq .action -r)
	time=$(echo "$line" | jq .time -r)
	data=$(echo "$line" | jq .data -r)

	# We'll track how long this action takes.
	actionStart=$(date +%s.%N)
	if [[ "$action" == "checkout" ]]; then
		git checkout -q "$data"
		echo "$(tput bold)Next step: $(git show --pretty=%s | head -n 1)$(tput sgr0)"
		sleep 2
		echo "File changed: $(git show --name-only | tail -n 1)"
		sleep 2
		git show --pretty | tail -n+9
		sleep 5
	elif [[ "$action" == "cmd" ]]; then
		pe "$data"
	else
		echo "???? ERROR"
		exit 42
	fi

	# Delay
	actionEnd=$(date +%s.%N)
	# This took
	actionObsTime=$(echo "$actionEnd - $actionStart" | bc -l)
	actionExpTime=$time
	# If this took less time than expected,

	res=$(echo "$actionObsTime < $actionExpTime" |bc -l);
	if (( res )); then
		toWait=$(echo "$actionExpTime - $actionObsTime" | bc -l)
		# Wait to catch up to that value
		sleep "$toWait"
	fi
	currentTime=$(echo "$actionEnd - $scriptStart" | bc -l)
	echo -e "$currentTime\t$action\t$data\t$time" >> "$SCRIPT_LOG"
done
IFS="$OLDIFS"
