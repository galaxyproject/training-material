#!/bin/bash
# Check the toolshed-y ones.
grep -n '\{% tool .*' -R topics | \
	python $(dirname $0)/check-broken-tool-links.py
