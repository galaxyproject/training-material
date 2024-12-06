#!/usr/bin/env bash
set -e

shed-tools install -g https://"${GALAXY_HOSTNAME}" -a adminkey --name minimap2 --owner iuc --section_label Mapping
galaxy-tool-test --galaxy-url https://"${GALAXY_HOSTNAME}" --key adminkey --tool-id minimap2 --test-index 0
