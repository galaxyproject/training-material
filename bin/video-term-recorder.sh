#!/bin/bash
OUR_DIR=$(dirname "$0")
SCRIPT="$1";
SCRIPT_LOG="$2";
ANSIBLE_HOST_OVERRIDE="$3";
GIT_GAT="$4";
rm -f "$SCRIPT_LOG"
# Disable git's paging which can cause a hang
export GIT_PAGER=cat
# Prevent pip from shouting everywhere
pip config --user set global.progress_bar off
# Setup the demo-magic
. ${OUR_DIR}/video-term-demo-magic.sh -n
# Go into right dir
cd ${GIT_GAT} || exit
# Hide our demo-magic activation
clear

OLDIFS="$IFS"
IFS=$'\n'
scriptStart=$(date +%s.%N)
for line in $(cat $SCRIPT | jq '.[]' -c); do
	# These aren't fast.
	action=$(echo "$line" | jq .action -r)
	time=$(echo "$line" | jq .time -r)
	data=$(echo "$line" | jq .data -r)

	# We'll track how long this action takes. This is also where we should start talking.
	actionStart=$(date +%s.%N)
	currentTime=$(echo "$actionStart - $scriptStart" | bc -l)
	echo -e "$currentTime\t$action" >> "$SCRIPT_LOG"
	if [[ "$action" == "checkout" ]]; then
		git checkout -q "$data"
		echo "$(tput bold)Next step: $(git show --pretty=%s | head -n 1)$(tput sgr0)"
		edited_file_name=$(git show --name-only | tail -n 1)
		sleep 1
		# Fake the PS1 to show the command they should use
		p "nano $edited_file_name"
		sleep 2
		# Checkout the previous commit so we can show the diff properly
		git checkout -q "$data^1"
		# This will pretend to edit it in """"vim""""
		python3 ${OUR_DIR}/video-diffplayer.py \
			--diff <(git show "$data") \
			--nosave \
			--speed 0.1 \
			--session-min-length "$time"
		# But now we actually need to be on the commit we say we are
		git checkout -q "$data"
	elif [[ "$action" == "cmd" ]]; then
		sed -i "s/gat.*.galaxyproject.eu/${ANSIBLE_HOST_OVERRIDE}/" ${GIT_GAT}/hosts
		echo 'password' > ${GIT_GAT}/.vault-password.txt
		pe "$data"
		git checkout -q -- hosts 2>/dev/null
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
done
IFS="$OLDIFS"
