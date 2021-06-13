#!/usr/bin/env bash
set -e

GALAXY_URL=https://"${GALAXY_HOSTNAME}" GALAXY_API_KEY=adminkey N_CHARS=14 N_THREADS=1 python3 $(dirname $0)/run_test.py
GALAXY_URL=https://"${GALAXY_HOSTNAME}" GALAXY_API_KEY=adminkey N_CHARS=17 N_THREADS=2 python3 $(dirname $0)/run_test.py
