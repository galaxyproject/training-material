# Install an entire topic to a Galaxy instance
#
# usages: install_topic_requirements.sh <path to topic> <galaxy url> <api key>

topic=$1
galaxy_url=$2
api_key=$3

for tutorial in ${topic}/tutorials/*
do
    echo "Installing Tutorial ${tutorial}"
    `dirname $0`/install_tutorial_requirements.sh $tutorial $galaxy_url $api_key
done

