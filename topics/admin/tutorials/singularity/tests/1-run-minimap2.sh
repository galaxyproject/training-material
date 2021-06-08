#!/usr/bin/env bash
set -e

shed-tools install -g https://"${GALAXY_HOSTNAME}" -a "$GALAXY_API_KEY" --name minimap2 --owner iuc --section_label Mapping
galaxy-tool-test --galaxy-url https://"${GALAXY_HOSTNAME}" --key "$GALAXY_API_KEY" --tool-id minimap2 --test-index 0
