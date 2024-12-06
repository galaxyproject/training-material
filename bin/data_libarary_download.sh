#!/bin/bash
set -e

# set API key and galaxy instance

API_KEY=${GALAXY_DEFAULT_ADMIN_KEY:-fakekey}
galaxy_instance="http://localhost:8080"

# setting up conda 
GALAXY_CONDA_PREFIX=/tool_deps/_conda export PATH=$GALAXY_CONDA_PREFIX/bin/:$PATH

if [ "$1" = '-d' ];
then
    if [ "$(head -n 1 /tutorials/data-library_all.yaml)" != "{}" ];
    then
        echo "Starting Galaxy in the background"
        #this is needed to use ephermeris
        startup &
        galaxy-wait -g $galaxy_instance
        echo "Downloading data from data-library_all.yaml"
        setup-data-libraries -i '/tutorials/data-library_all.yaml' -g $galaxy_instance -a $API_KEY -v
        #because galaxy runs in the background, it stops running after the data is downloaded. A second startupcommand is needed:
        startup
    else
        echo "No data libraries to install"
        echo "Starting Galaxy"
        startup
    fi

else
    echo "Use '-d' parameter to download tutorial data."
    echo "Starting Galaxy"
    startup
fi