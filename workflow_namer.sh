#!/bin/bash
for topicdir in ./topics/*
do
    topic=$(basename $topicdir)
    echo "TOPIC: ${topic^^} "
    echo "----------------------------------------------" 
    for tutdir in $topicdir/tutorials/*
    do
        tut=$(basename $tutdir)
        # install tools and workflows
        if [ -d $tutdir/workflows/ ];
        then
            for w in $tutdir/workflows/*.ga
            do
                python3 ./workflow_namer.py $w $topic
            done
        fi
    done
done