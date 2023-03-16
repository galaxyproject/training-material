#!/bin/bash
DIR=/tmp/git-gat
op="$1"

declare -a tutorials
tutorials=(admin/ansible-galaxy admin/tus admin/cvmfs admin/singularity admin/tool-management admin/data-library admin/connect-to-compute-cluster admin/job-destinations admin/pulsar admin/gxadmin admin/monitoring admin/tiaas admin/reports admin/ftp)
#tutorials=(admin/wireguard-headscale)
#tutorials=(admin/wireguard)

#echo "${tutorials[0]}"
#exit 1;
if [[ "$op" == "export" ]]; then
	rm -rf "${DIR}"
	mkdir -p "${DIR}"
	mkdir -p "${DIR}/.scripts"

	# Setup readme as the root commit

	cat > ${DIR}/0-commit-0000-root-commit.patch <<-EOF
		From: The Galaxy Training Network <galaxytrainingnetwork@gmail.com>
		Date: Mon, 15 Feb 2021 14:06:56 +0100
		Subject: admin/init/0000: Root commit


		--- /dev/null
		+++ b/readme.md
		@@ -0,0 +1,11 @@
		+# GIT-GAT
		+
		+This is a git repository with the current <abbr title="Galaxy Admin Training">GAT</abbr> history. See the current [GAT schedule](https://gxy.io/gat).
		+
		+This is built from [the GTN's library of admin training](https://training.galaxyproject.org/training-material/topics/admin/)
		+
		+Extra | Data
		+--- | ---
		+Date | $(date --rfc-3339=seconds)
		+Github Run ID | [$GITHUB_RUN_ID](https://github.com/galaxyproject/training-material/actions/runs/$GITHUB_RUN_ID)
		+GTN Commit | [$GITHUB_SHA](https://github.com/galaxyproject/training-material/tree/$GITHUB_SHA)
		--
		2.25.1
	EOF


	# Then do the rest
	for idx in "${!tutorials[@]}"; do
		echo "Processing ${tutorials[$idx]}"
		folder=$(echo "${tutorials[$idx]}" | cut -d / -f 1)
		tuto=$(echo "${tutorials[$idx]}" | cut -d / -f 2)

		python3 bin/knit-frog.py \
			topics/${folder}/tutorials/${tuto}/tutorial.md \
			${DIR}/$(( idx + 10 ))-${tuto};
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
	# Export everything BUT the root commit
	git format-patch "${root_commit}..HEAD"

	# Go back
	cd "${CURRENT_DIR}" || exit

	# Import all of the patches
	for idx in "${!tutorials[@]}"; do
		if [[ "${tutorials[$idx]}" != "admin/tool-management" ]]; then
		folder=$(echo "${tutorials[$idx]}" | cut -d / -f 1)
		tuto=$(echo "${tutorials[$idx]}" | cut -d / -f 2)

		python3 bin/knit.py \
			topics/${folder}/tutorials/${tuto}/tutorial.md \
			--patches ${DIR}/*${folder}-${tuto}*.patch
		fi
	done
elif [[ "$op" == "deploy" ]]; then
	if [[ ! -d "${DIR}" ]]; then
		echo "Error, ${DIR} is missing"
		exit 1
	fi

	cd ${DIR} || exit
	git init && \
		git add .scripts/
		git commit -a -m 'Add scripts directory'
		git am -3 -- *.patch && \
		git remote add origin git@github.com:hexylena/git-gat.git && \
		git push -f origin
elif [[ "$op" == "roundtrip" ]]; then
	rm -rf ${DIR}
	bash $0 export
	cd ${DIR} || exit
	git init && \
		git am -3 -C2 -- *.patch
	cd -
	bash $0 import
else
	echo "$0 <import|export|deploy|roundtrip>"
	exit 1;
fi


