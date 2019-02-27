#! /bin/bash

# Script to install a tutorial to a running Galaxy instance
#
# usage: install_tutorial.sh </path/to/tutorial> <galaxy_url> <admin_api_keyoexample.org
#
# example: sh bin/install_tutorial.sh topics/metagenomics/tutorials/mothur-miseq-sop http://localhost:8080 admin
#
# make sure you have ephemeris installed:
#    pip install ephemeris -U
#

tutorial=$1
galaxy_url=$2
api_key=$3

    
# install workflows and tools
echo "Installing workflows"
if [ -d ${tutorial}/workflows ]
then
    echo " - Extracting tools from workflows"
    for w in ${tutorial}/workflows/*.ga
    do
        workflow-to-tools -w $w -o ${tutorial}/workflows/wftools.yaml -l $(basename ${tutorial})
        echo " - Installing tools from workflow $(basename $w)" 
        shed-tools install -g ${galaxy_url} -a ${api_key} -t ${tutorial}/workflows/wftools.yaml
        rm ${tutorial}/workflows/wftools.yaml
    done
    echo " - Installing workflows"
    workflow-install --publish_workflows -g ${galaxy_url} -a ${api_key} -w ${tutorial}/workflows/

else
    echo "No workflows to install (no directory named workflow present)"
fi

# install data libraries
echo "Populating Data Libraries"
if [ -f ${tutorial}/data-library.yaml ]
then
    setup-data-libraries -g ${galaxy_url} -a ${api_key} -i ${tutorial}/data-library.yaml
else
    echo "No data library to install (no file named data-libraries.yaml present)"
fi

# install reference data
echo "Installing reference data"
if [ -f ${tutorial}/data-manager.yaml ]
then
    run-data-managers -g ${galaxy_url} -a ${api_key} --config ${tutorial}/data-manager.yaml

else
    echo "No reference data to install (no file named data-manager.yaml present)"
fi

# install tours --> not possible via API yet
