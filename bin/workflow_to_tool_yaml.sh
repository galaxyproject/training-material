#!/bin/bash

# extract the tools.yaml file from the workflows

for topic in topics/*
do
    echo "Topic: $topic"
    for tutorial in "$topic/tutorials/"*
    do
        echo " - Tutorial: $tutorial"
        if [ -d "${tutorial}/workflows" ]
        then
            wf_files=$(find "${tutorial}/workflows/" -name '*.ga' | paste -sd " ")
            num_wf_files=$(find "${tutorial}/workflows/" -name '*.ga' | wc -l)
            if (( num_wf_files > 0 )); then
                workflow-to-tools -w ${wf_files} -o "${tutorial}/tools.yaml" -l "$(basename "${tutorial}")"
                echo "    ..done. Tool list in ${tutorial}/tools.yaml"
            else
                echo "    No workflows to install (no files found)"
            fi
        else
            echo "    No workflows to install (no directory named workflow present)"
        fi
    done
done
