#!/bin/bash
DIR=/tmp/gat-git
op="$1"
if [[ "$op" == "export" ]]; then
	mkdir -p ${DIR}

	idx=0
	for i in {ansible-galaxy,singularity,cvmfs,data-library,connect-to-compute-cluster,job-destinations,pulsar}; do
		python bin/knit-frog.py \
			topics/admin/tutorials/$i/tutorial.md \
			${DIR}/${idx}-$i;

		idx=$(( idx + 1 ));
	done
elif [[ "$op" == "import" ]]; then
	if [[ ! -d "${DIR}" ]]; then
		echo "Error, ${DIR} is missing"
		exit 1
	fi

	# Store current dir
	CURRENT_DIR=$(pwd)
	cd "${DIR}" || exit

	# Export root commit
	root_commit=$(git log --pretty=oneline | tail -n 1 | cut -f1 -d' ')
	git format-patch --root "$root_commit"
	# Export the rest
	git format-patch "${root_commit}..HEAD"

	# Go back
	cd "${CURRENT_DIR}" || exit

	# Import all of the patches
	for i in {ansible-galaxy,singularity,cvmfs,data-library,connect-to-compute-cluster,job-destinations,pulsar}; do
		python bin/knit.py \
			topics/admin/tutorials/$i/tutorial.md \
			--patches ${DIR}/*admin-${i}*.patch
	done
else
	echo "$0 <import|export>"
	exit 1;
fi


