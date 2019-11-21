#!/bin/bash
set -e

if [ "$1" = '-d' ];
then
    if [ "$(head -n 1 /tutorials/data-library_all.yaml)" != "{}" ];
    then
        echo "Starting Galaxy in the background"
        #this is needed to use ephermeris
        startup &
        galaxy_instance="http://localhost:8080"
        galaxy-wait -g $galaxy_instance
        echo "Downloading data from data-library_all.yaml"
        setup-data-libraries -i '/tutorials/data-library_all.yaml' -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD -v
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