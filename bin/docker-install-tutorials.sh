# install everything for all tutorials of a topic

set -e
galaxy_instance="http://localhost:8080"

# launch the instance
echo " - Launch a lite instance"
startup_lite
sleep 120 #galaxy-wait --timeout 60

# install other tutorial materials
for dir in /tutorials/*
do
    echo "For tutorial: $dir"

    # install tools
    echo " - Installing tools"
    shed-install -t $dir/tools.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD

    # install workflows (TODO: make them shared instead of just under admin user account?)
    if [ -d $dir/workflows/ ];
    then
        echo " - Installing workflows"
        workflow-install --workflow_path $dir/workflows/ -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo "No workflow to install"
    fi

    # install data libraries
    echo " - Installing data libraries"
    setup-data-libraries -i $dir/data-library.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD

    # install reference data? (discussion: do this at build or run time?)
    if [ -f $dir/data-manager.yaml ]
    then
        echo " - Installing reference data"
        run-data-managers --config $dir/data-manager.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo "No reference data to install (no data-manager.yaml file present)"
    fi

    # install tours
    dir_name="$(dirname $dir)"
    dir_name="$(basename $dir_name)"
    for t in $dir/tours/*
    do
        # prefix tour file name with tutorial name to avoid clashes
        fname=$dir_name-$(basename $t)
        echo " - Installing tour: $t as $fname"
        cp $t $GALAXY_ROOT/config/plugins/tours/$fname
    done
done
