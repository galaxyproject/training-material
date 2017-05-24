# Script to install a tutorial to a running Galaxy instance
#
# usage: install_tutorial.sh <topic_name> <galaxy_url> <admin_api_key>
#
# example: sh bin/install_tutorial.sh metagenomics http://localhost:8080 admin
#
# make sure you have ephemeris installed:
#    pip install ephemeris -U
#

topic=$1
galaxy_url=$2
api_key=$3

# install tools
shed-install -g ${galaxy_url} -a ${api_key} -t ${topic}/tools.yaml

# install data libraries
setup-data-libraries -g ${galaxy_url} -a ${api_key} -i ${topic}/data_libraries.yaml

# install workflows
workflow-install -g ${galaxy_url} -a ${api_key} -w ${topic}/workflows/

# run data managers
run-data-managers -g ${galaxy_url} -a ${api_key} --config ${topic}/data_manager.yaml

# install tours --> not possible via API yet
