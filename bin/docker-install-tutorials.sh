# install everything for all tutorials of a topic

galaxy_instance="http://localhost:8080"

# start galaxy
startup_lite
/galaxy-sleep.py $galaxy_instance # replace with ephemeris sleep once updated in pip

# install tutorial materials
cd /tutorials
for dir in *
do
    echo "Installing tutorial: $dir"

    # install tools
    echo " - Installing tools"
    shed-install -t $dir/tools.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD

    # install workflows (TODO: make them shared instead of just under admin user account?)
    echo " - Installing workflows"
    workflow-install --workflow_path $dir/workflows/ -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD

    # install data libraries
    echo " - Installing data libraries"
    setup-data-libraries -i $dir/data-library.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD

    # install tours
    for t in $dir/tours/*
    do
        # prefix tour file name with tutorial name to avoid clashes
        fname=$dir-$(basename $t)
        echo " - Installing tour: $t as $fname"
        cp $t $GALAXY_ROOT/config/plugins/tours/$fname
    done

    # install reference data? (discussion: do this at build or run time?)
    if [ -f $dir/data-manager.yaml ]
    then
        echo "Installing reference data"
        run-data-managers --config $dir/data-manager.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo "No reference data to install (no data-manager.yaml file present)"
    fi
done
