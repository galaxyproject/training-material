# install everything for all tutorials of a topic

set -e
galaxy_instance="http://localhost:8080"

# launch the instance
echo " - Starting Galaxy.. \n"
startup_lite

# wait until galaxy has started
galaxy-wait -g $galaxy_instance

# install other tutorial materials
for dir in /tutorials/*
do
    echo "Installing tutorial: $dir"

    # install tools
    if [ -f $dir/tools.yaml ]
    then
        echo " - Installing tools"
        shed-install -t $dir/tools.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo " - No tools to install (no file named tools.yaml present)"
    fi

    # install workflows (TODO: make them shared instead of just under admin user account?)
    if [ -d $dir/workflows/ ];
    then
        echo " - Installing workflows"
        workflow-install --workflow_path $dir/workflows/ -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo " - No workflows to install (no directory named workflows present)"
    fi

    # install data libraries
    if [ -f $dir/data-library.yaml ]
    then
        echo " - Installing data libraries"
        setup-data-libraries -i $dir/data-library.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo " - No data libraries to install (no file named data-library.yaml present)"
    fi

    # install reference data? (discussion: do this at build or run time?)
    if [ -f $dir/data-manager.yaml ]
    then
        echo " - Installing reference data"
        run-data-managers --config $dir/data-manager.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
    else
        echo " - No reference data to install (no file named data-manager.yaml present)"
    fi

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
