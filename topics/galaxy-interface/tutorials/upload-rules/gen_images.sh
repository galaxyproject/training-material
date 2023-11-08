#!/bin/bash
# Usage:
#    bash topics/introduction/tutorials/galaxy-intro-rules/gen_images.sh /path/to/galaxy
#
# Description:
# This is meant to show how the images were generated more than for end user consumption
# so things are a bit messy and sub-optimal.

export CWD=`pwd`
export THIS_DIR="$(dirname "$0")"

GALAXY_ROOT="$1"
cd $GALAXY_ROOT

# in case THIS_DIR is relative go back to where we started
make client
GALAXY_SKIP_CLIENT_BUILD=1 GALAXY_RUN_WITH_TEST_TOOLS=1 GALAXY_CONFIG_OVERRIDE_DATABASE_CONNECTION=postgres://postgres:mysecretpassword@localhost:5432/galaxy sh run.sh &

sleep 60

# Set Galaxy Selenium tests to generate screenshots into a clean directory.
rm -rf screens
export GALAXY_TEST_SCREENSHOTS_DIRECTORY=`pwd`/screens

. .venv/bin/activate

# Run Galaxy Selenium tests corresponding to the examples.
GALAXY_TEST_EXTERNAL=http://localhost:8080/ nosetests test/selenium_tests/test_uploads.py:UploadsTestCase.test_rules_example_1_datasets
GALAXY_TEST_EXTERNAL=http://localhost:8080/ nosetests test/selenium_tests/test_uploads.py:UploadsTestCase.test_rules_example_2_list
GALAXY_TEST_EXTERNAL=http://localhost:8080/ nosetests test/selenium_tests/test_uploads.py:UploadsTestCase.test_rules_example_3_list_pairs
GALAXY_TEST_EXTERNAL=http://localhost:8080/ nosetests test/selenium_tests/test_uploads.py:UploadsTestCase.test_rules_example_4_accessions
GALAXY_TEST_EXTERNAL=http://localhost:8080/ nosetests test/selenium_tests/test_uploads.py:UploadsTestCase.test_rules_example_5_matching_collections
GALAXY_TEST_EXTERNAL=http://localhost:8080/ nosetests test/selenium_tests/test_uploads.py:UploadsTestCase.test_rules_example_6_nested_lists

cd "$CWD"
cd "$THIS_DIR/../.."
rm -f images/rules/rules_example*
cp "$GALAXY_ROOT"/screens/rules_example* images/rules
