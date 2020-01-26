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
for tutdir in $topicdir/tutorials/*
do
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
                shed-tools install -t $tutdir/workflows/wftools.yaml -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD -v && break
                n=$[$n+1]
                sleep 5
                echo " - Retrying shed-tools install "
            done        
            rm $tutdir/workflows/wftools.yaml          
        done
        echo " - Installing workflows"
        workflow-install --publish_workflows --workflow_path $tutdir/workflows/ -g $galaxy_instance -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD
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

    echo "Finished installation of $tut tutorial \n"
done

cd /tutorials/
python /mergeyaml.py > ./data-library_all.yaml