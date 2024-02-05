#!/bin/bash
OUR_DIR=$(dirname $(realpath $0))
SPEED=0.01
# Disable git's paging which can cause a hang
export GIT_PAGER=cat

# Prevent pip from shouting everywhere
pip config --user set global.progress_bar off
# Setup the demo-magic
. ${OUR_DIR}/video-term-demo-magic.sh -n
# Hide our demo-magic activation
${OUR_DIR}/knit-automated-dev.sh export
# Setup our git dir
cd /tmp/tools-from-scratch 
git init 
git config --local --add commit.gpgsign false
git am -3 *
# Clear the screen
clear

for commit in $(git log --pretty=oneline | head -n-1 | tac | cut -f 1 -d' '); do
	git checkout -q $commit^
	echo "$(tput bold)Next step: $(git show --pretty=%s $commit | head -n 1)$(tput sgr0)"
	edited_file_name=$(git show --name-only $commit| tail -n 1)
	sleep 1
	# Fake the PS1 to show the command they should use
	p "nano $edited_file_name"
	sleep 2
	# Checkout the previous commit so we can show the diff properly
	# This will pretend to edit it in """"vim""""
	python3 ${OUR_DIR}/video-diffplayer.py \
		--diff <(git show "$commit") \
		--nosave \
		--speed ${SPEED}
done
