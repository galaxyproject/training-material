#!/bin/bash

# extract the tools.yaml file from the workflows

for topic in topics/*
do
    echo "Topic: $topic"
    for tutorial in $topic/tutorials/*
    do
        echo " - Tutorial: $tutorial"
        if [ -d ${tutorial}/workflows ]
        then
            workflow-to-tools -w ${tutorial}/workflows/*.ga -o ${tutorial}/tools.yaml -l $(basename ${tutorial})
            echo "    ..done. Tool list in ${tutorial}/tools.yaml"
        else
            echo "    No workflows to install (no directory named workflow present)"
        fi
    done
done
