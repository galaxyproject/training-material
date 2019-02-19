# install everything for all tutorials of a topic

set -e
galaxy_instance="http://localhost:8080"

# launch the instance
echo " - Starting Galaxy.. \n"

export GALAXY_CONFIG_TOOL_PATH=/galaxy-central/tools/
startup_lite

# wait until galaxy has started
galaxy-wait -g $galaxy_instance

# run tutorial install as user galaxy
su - $GALAXY_USER

# install other tutorial materials
for dir in /tutorials/*
do
    echo "Installing tutorial: $dir"

    # install tools and workflows
    if [ -d $dir/workflows/ ];
    then
        echo " - Extracting tools from workflows"
        for w in $dir/workflows/*.ga
        do
            workflow-to-tools -w $w -o $dir/workflows/wftools.yaml -l $(basename $dir)
            echo " - Installing tools from workflow $(basename $w)" 
            shed-tools install -t $dir/workflows/wftools.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
            rm $dir/workflows/wftools.yaml
        done
        echo " - Installing workflows"
        workflow-install --publish_workflows --workflow_path $dir/workflows/ -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo " - No workflows to install (no directory named workflows present)"
    fi

    # install reference data? (discussion: do this at build or run time?)
    # We are using CVMFS for the moment.
    #if [ -f $dir/data-manager.yaml ]
    #then
    #    echo " - Installing reference data"
    #    run-data-managers --config $dir/data-manager.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    #else
    #    echo " - No reference data to install (no file named data-manager.yaml present)"
    #fi

    # install tours
    dir_name="$(dirname $dir)"
    dir_name="$(basename $dir_name)"
    if [ -d $dir/tours/ ];
    then
        echo " - Installing tours"
        for t in $dir/tours/*
        do
            # prefix tour file name with tutorial name to avoid clashes
            fname=$dir_name-$(basename $t)
            echo "   - Installing tour: $t as $fname"
            cp $t $GALAXY_ROOT/config/plugins/tours/$fname
        done
    else
        echo " - No tours to install (no directory named tours present)"
    fi

    echo "Finished installation of $dir tutorial \n"
done

cd /tutorials/
python /mergeyaml.py > ./data-library_all.yaml
# check if data-library_all.yaml is not empty
if [ "$(head -n 1 ./data-library_all.yaml)" != "{}" ];
then
    setup-data-libraries -i ./data-library_all.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD -v
else
    echo "No data libraries to install"
fi
