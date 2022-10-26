#!/bin/bash

# set API key
API_KEY=${GALAXY_DEFAULT_ADMIN_KEY:-fakekey}

# Enable Test Tool Shed
echo ".. Setting up conda"
export GALAXY_CONFIG_TOOL_SHEDS_CONFIG_FILE=$GALAXY_HOME/tool_sheds_conf.xml

. /tool_deps/_conda/etc/profile.d/conda.sh
conda activate base

if pgrep "supervisord" > /dev/null
then
    echo ".. System is up and running. Starting with the installation."
    export PORT=80
else
    # start Database
    echo ".. Starting database"
    export PORT=8080
    service postgresql start
    install_log='galaxy_install.log'

    # wait for database to finish starting up
    STATUS=$(psql 2>&1)
    while [[ ${STATUS} =~ "starting up" ]]
    do
      echo ".. Waiting for database: $STATUS"
      STATUS=$(psql 2>&1)
      sleep 1
    done
    # start Galaxy
    echo ".. Starting Galaxy"
    # Unset SUDO_* vars otherwise conda run chown based on that
    sudo -E -u galaxy -- bash -c "unset SUDO_UID; \
        unset SUDO_GID; \
        unset SUDO_COMMAND; \
        unset SUDO_USER; \
        ./run.sh -d $install_log --pidfile galaxy_install.pid --http-timeout 3000"

    galaxy_install_pid=`cat galaxy_install.pid`
    galaxy-wait -g http://localhost:$PORT -v --timeout 120
    echo ".. Galaxy is running"
fi

# Create the admin user if not already done
# Starting with 20.05 this user is only created at first startup of galaxy
# We need to create it here for Galaxy Flavors = installing from Dockerfile
if [[ ! -z $GALAXY_DEFAULT_ADMIN_USER ]]
    then
        (
        cd $GALAXY_ROOT
        . $GALAXY_VIRTUAL_ENV/bin/activate
        echo ".. Creating admin user $GALAXY_DEFAULT_ADMIN_USER with key $GALAXY_DEFAULT_ADMIN_KEY and password $GALAXY_DEFAULT_ADMIN_PASSWORD if not existing"
        python /usr/local/bin/create_galaxy_user.py --user "$GALAXY_DEFAULT_ADMIN_EMAIL" --password "$GALAXY_DEFAULT_ADMIN_PASSWORD" \
        -c "$GALAXY_CONFIG_FILE" --username "$GALAXY_DEFAULT_ADMIN_USER" --key "$GALAXY_DEFAULT_ADMIN_KEY"
        )
fi

# install tutorial materials
echo ".. Starting installation of the tutorials."
for tutdir in $topicdir/tutorials/*
do
    echo "-------------------------------------------------------------"
    tut=$(basename $tutdir)
    echo "Installing tutorial: $tut"
    # install tools and workflows
    if [ -d $tutdir/workflows/ ];
    then
        echo " - Extracting tools from workflows"
        for w in $tutdir/workflows/*.ga
        do
            workflow-to-tools -w $w -o $tutdir/workflows/wftools.yaml -l $tut
            echo " - Installing tools from workflow $(basename $w)" 
            n=0
            until [ $n -ge 3 ]
            do
                shed-tools install -t $tutdir/workflows/wftools.yaml -g "http://localhost:$PORT" -a $API_KEY && break
                n=$[$n+1]
                sleep 5
                echo " - Retrying shed-tools install "
            done        
            rm $tutdir/workflows/wftools.yaml          
        done
        echo " - Installing workflows"
        workflow-install --publish_workflows --workflow_path $tutdir/workflows/ -g "http://localhost:$PORT" -a $API_KEY
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
    if [ -d $tutdir/tours/ ];
    then
        echo " - Installing tours"
        for t in $tutdir/tours/*
        do
            fname=$tut-$(basename $t)
            echo "   - Installing tour: $t as $fname"
            cp $t $GALAXY_ROOT/config/plugins/tours/$fname
        done
    else
        echo " - No tours to install (no directory named tours present)"
    fi

    echo " - Finished installation of $tut tutorial"
done

echo "-------------------------------------------------------------"
echo ".. Generating data-library_all.yaml file"
cd /tutorials/
python /mergeyaml.py > ./data-library_all.yaml

exit_code=$?

if [ $exit_code != 0 ] ; then
    if [ "$2" == "-v" ] ; then
        echo ".. Installation failed, Galaxy server log:"
        cat $install_log
    fi
    exit $exit_code
fi

if ! pgrep "supervisord" > /dev/null
then
    echo ".. Shutting down Galaxy and postgresql"
    # stop everything
    sudo -E -u galaxy /galaxy-central/run.sh --stop --pidfile /galaxy-central/galaxy_install.pid
    rm /galaxy-central/$install_log
    service postgresql stop
fi

echo ".. Installation is finished"
