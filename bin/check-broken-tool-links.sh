#!/bin/bash
set -e -o pipefail
# Check the toolshed-y ones.
grep -n '\{% tool .*' -R topics | \
	python3 $(dirname $0)/check-broken-tool-links.py
