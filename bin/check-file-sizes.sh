#!/bin/bash

echo "checking file sizes"
exitcode=0
# check that test_data folder for workflows do not exceed 10MB

for TEST_DATA in topics/*/tutorials/*/workflows/test-data; do

if [[ $(du -sb $TEST_DATA | cut -f1) -gt 10485760 ]]; then
  if [[ $exitcode -eq 0 ]]; then
    echo "ERROR: test-data folder may not exceed 10MB."
    echo "  To reduce size, please use zenodo links for files (and replace 'path:' with 'location:'). And consider replacing output files with asserts"
  fi
  echo "  Folder too big: $(du -sh $TEST_DATA)"
  exitcode=1
fi

done

exit $exitcode
