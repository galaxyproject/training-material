#!/bin/bash

# extract the tools.yaml file from the workflows

for topic in topics/*
do
    echo "Topic: $topic"
    for tutorial in $topic/tutorials/*
    do
        echo "Tutorial: $tutorial"
        if [ -d ${tutorial}/workflows ]
        then
            for w in ${tutorial}/workflows/*.ga
            do
                workflow-to-tools -w $w -o ${tutorial}/tools.yaml -l $(basename ${tutorial})
            done
        else
            echo "No workflows to install (no directory named workflow present)"
        fi
    done
done
