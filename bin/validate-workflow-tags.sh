#!/bin/bash
exit_with=0

virtualenv gxformat2
. gxformat2/bin/activate
pip install gxformat2

for topicdir in ./topics/*
do
    topic=$(basename $topicdir)
    for tutdir in $topicdir/tutorials/*
    do
        if [ -d $tutdir/workflows/ ];
        then
            for w in $tutdir/workflows/*.ga
            do
                echo "Checking $w"
                if gxwf-lint --training-topic $topic $w;
                then
                    echo "-------------Invalid workflow--------------"
		            exit_with=1
                fi
            done
        fi
    done
done

exit $exit_with
