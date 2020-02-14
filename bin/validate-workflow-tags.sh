#!/bin/bash
exit_with=0

for topicdir in ./topics/*
do
    topic=$(basename $topicdir)
    for tutdir in $topicdir/tutorials/*
    do
        if [ -d $tutdir/workflows/ ];
        then
            for w in $tutdir/workflows/*.ga
            do
                if [ -f $w ]
                then
                    echo "Checking $w"
                    if planemo workflow_lint $w --topic_name $topic;
                    then
                        exit_with=1
                    fi
                fi
            done
        fi
    done
done

exit $exit_with
